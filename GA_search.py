import pygad
import subprocess
import pandas as pd
import os
import numpy as np
import nibabel as nib
import glob
import json
import time
import random
import string
from pathlib import Path
import csv
import math
import scipy.io as scio


# value boundaries
# fa_threshold [0,0.7]
# min_length [0.1, 5]
# max_len [10,200]
# turning_angle [15, 60]
# step_size [0.01, 0.05] step size is in mm. our data has a resolution of 0.025 mm. ranges from sub-voxel to 2-voxels
# otsu_threshold [0.3, 1.2]
DEBUG = True
CLUSTER = True
SOL_PER_GENERATION = 24
if DEBUG:
    SOL_PER_GENERATION = 3
    SOL_PER_GENERATION = 10

class GA_pipeline:
    def __init__(self, ngen, experiment_table_path, runno_list, debug, dsi_studio, exclusion_list, num_parents_mating, gene_space, gene_type, project_code):
        # HARRISON ADDED RANDOMLY CHOSEN VALUE (previously unset)
        self.project_code = project_code
        self.debug = debug
        self.exclusion_list = exclusion_list
        self.runno_list = runno_list
        self.experiment_table_path = experiment_table_path
        self.ngen = ngen
        self.num_parents_mating = num_parents_mating
        self.num_genes = 5
        self.gene_space = gene_space
        self.gene_type = gene_type
        self.best_sol_uids = []
      
        # UPDATED VERSION
        self.output_dir_base = "{}/{}/{}".format(os.environ["BIGGUS_DISKUS"], self.project_code, "genetic_parameter_sets")
        if DEBUG:
            self.output_dir_base = "{}/{}/debug_test_08-11-2024/{}".format(os.environ["BIGGUS_DISKUS"], self.project_code, "genetic_parameter_sets")
            self.output_dir_base = "{}/{}/debug_test_10-by-25/{}".format(os.environ["BIGGUS_DISKUS"], self.project_code, "genetic_parameter_sets")
        self.src_sir = "{}/../src".format(self.output_dir_base)
        self.fib_sir = "{}/../fib".format(self.output_dir_base)
        # this one is only used in the cluster. otherwise, pull from archive
        # TODO: it seems like the cluster archive connection is working again, test and unify solution (if so)
        self.label_dir = "{}/../labels".format(self.output_dir_base)
        self.bash_wrapper_dsi_studio_dir = "{}/../bash_wrapper_dsi_studio".format(self.output_dir_base)
        self.bash_wrapper_omni_manova_dir = "{}/../bash_wrapper_omni_manova".format(self.output_dir_base)
        self.omni_manova_dir_base = "{}/../genetic_search".format(self.output_dir_base)
       
        # TODO: is there a way to get this archive base path from the environment?
        self.archive_dir = "A:/{}/research".format(self.project_code)
        if CLUSTER:
            self.archive_dir = "/mnt/nclin-comp-pri.dhe.duke.edu/dusom_civm-atlas/{}/research/".format(self.project_code) 

        self.matlabb = "matlab"
        self.dsi_studio = dsi_studio


    def setup_dsi_studio_trk_call(self, experiment: dict, fib_file, qa_nifti_file, label_file, output_dir):
        options_str = ""
        dsi_exclusion_list = ["uid", "fa_threshold", "output", "source", "connectivity", "connectivity_output"]
        for key in experiment:
            if key in dsi_exclusion_list:
                continue
            options_str += "--{}={} ".format(key, experiment[key])
        # inputs and outputs, these are different for every individual, so handled separately from the options dict
        options_str += "--output={} ".format(output_dir)
        options_str += "--source={} ".format(fib_file)
        options_str += "--connectivity={} ".format(label_file)
        # if we do not set this, then it will be saved right next to the trk file. we want this.
        # options_str += "--connectivity_output={} ".format(output_dir)

        threshold = self.find_vol_pct_threshold(qa_nifti_file, experiment["fa_threshold"])
        options_str += "--fa_threshold={} ".format(threshold)

        s = "{} --action=trk {}".format(self.dsi_studio, options_str)
        # print(s)
        return s


    def make_cluster_command(self, cmd, out_dir, job_name, memory="40G"):
        # this will convert a command to a cluster command
            # this is done either by srun --mem $cmd to immediately run a command
            # or sbatch --mem path/to/scrp.bash to schedule the command to run whenever
        # and then save that command to a bash file along with #!bash
        # returns the sbatch command
        #cmd = "srun --mem={} {}".format(memory, cmd)
        bash_wrapper = "{}/{}-{}.bash".format(out_dir, job_name, time.time())

        with open(bash_wrapper, "w") as f:
            f.write("#!/usr/bin/env bash\n")
            f.write(cmd)
        return "sbatch --mem={} {}".format(memory, bash_wrapper)




    def find_vol_pct_threshold(self, image, percent):
        nifti = nib.load(image)
        data, header = nifti.get_fdata(), nifti.header
        data = data[data != 0]
        data = np.sort(data)
        index = round((percent / 100) * len(data))
        threshold = data[index]
        print(threshold)
        return threshold


    def run_dsi_studio(self, experiment_list): # this changes
        '''csv_table = pd.read_csv(self.experiment_table_path, header=1, delimiter="\t")
        keys = csv_table.columns.tolist()
        vals = csv_table.iloc[0].tolist()
        const_keys = list(set(keys) - set(keys[1:13]))
        solution_keys = keys[1:13]
        const_vals = vals[0].append(vals[13:])

        experiment_list = dict(zip(const_keys, const_vals))
        sol_dic = dict(zip(solution_keys, sol))
        experiment_list.update(sol_dic)'''
        count = 0
        cmds = []
        for experiment in experiment_list:
            if self.debug:
                count += 1
            experiment_dir = "{}/{}".format(self.output_dir_base, experiment["uid"])
            if not os.path.exists(experiment_dir):
                os.makedirs(experiment_dir)
            process_list = []
            for runno in self.runno_list:
                output_dir = "{}/{}".format(experiment_dir, runno)
                if not os.path.exists(output_dir):
                    os.mkdir(output_dir)
                # i was having archive connection troubles on the cluster, so I just loaded in the labels file to a local directory
                # labels (i think) is all that was needed from the archive for this
                #label_file = "{}/connectome{}dsi_studio/labels/RCCF/{}_RCCF_labels.nii.gz".format(self.archive_dir, runno, runno)
                label_file = "{}/{}_RCCF_labels.nii.gz".format(self.label_dir, runno)
                fib_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz".format(self.fib_dir, runno)
                qa_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz.qa.nii.gz".format(self.fib_dir, runno)
                # check for completion by looking for any trk file in the outputs folder
                found_trk = glob.glob("{}/*.trk.gz".format(output_dir))
                found_tt = glob.glob("{}/*.tt.gz".format(output_dir))
                if len(found_trk) > 0 or len(found_tt) > 0:
                    print("work already complete for experiment {} runno {}".format(experiment["uid"], runno))
                    continue
                cmd = self.setup_dsi_studio_trk_call(experiment, fib_file, qa_file, label_file, output_dir)
                if CLUSTER:
                    # then we make a list of commands to be scheduled
                    # make a list and schedule them together because easier to be able to wait for all these jobs as a group
                    # that way we know when ready to go on to omni-manova
                    job_name = "connectome_{}_exp-{}".format(runno, experiment["uid"])
                    cmd = self.make_cluster_command(cmd, self.bash_wrapper_dsi_studio_dir, job_name)
                cmds.append(cmd)
            with open('{}/experiment.json'.format(experiment_dir), 'w') as fp:
                json.dump(experiment, fp)
            # i forget what these exit codes are for
            #exit_codes = [p.wait() for p in process_list]
        # now handle all the jobs
        if CLUSTER:
            print("sleeping 5 seconds to ensure all bash stub files are created")
            time.sleep(5)
            print("running all via sbatch")
            print(cmds)
            job_ids = []
            for cmd in cmds:
                proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE)
                # stdout is a buffered reader (bytes string)
                # must first decode it, and then split on spaces to get the last value (the job id)
                jid = proc.stdout.read().decode("utf-8").split(" ")[-1]
                # to strip off the trailing newline character
                jid = jid.strip("\n")
                print("sbatch job id: {}".format(jid))
                job_ids.append(jid)
            
            # call the pipeline_utilities 'cluster_wait' function, pass it list of slurm ids
            # turn my job ids list into a space-delimited string
            job_ids = " ".join(job_ids)
            cluster_wait_cmd = "pipeline_utilities cluster_wait {}".format(job_ids)
            print("waiting for dsi studio calls to complete")
            subprocess.run(cluster_wait_cmd.split(" "))
        else:
            for cmd in cmds:
                subprocess.run(cmd, shell=True)
        print("done scheduling DSI Studio calls")



    def run_omnimanova(self, gen):
        out_dir = "{}/Omni_Manova-{}".format(self.omni_manova_dir_base, gen)
        data_frame_path = "{}/dataframe.csv".format(out_dir)
        mat_script_template = "{}/analysis/Omni_Manova/Run_File/template_prototype_run_from_python.m".format(os.environ["WORKSTATION_CODE"])
        # MUST be turned OFF when running on the cluster. R doesn't work there. Optional to run on CITRIX
        run_R_analysis = 0;
        test_criteria = "{{'group1','group2','subgroup02','subgroup03', 'subgroup04', 'subgroup05','subgroup06','subgroup07','subgroup08','subgroup09','subgroup10','subgroup11','subgroup12'}}"

        timestamp = int(time.time())
        log_file = "{}/omni_manova_{}.log".format(out_dir, timestamp)
        # made a change to ALWAYS USE A UNIQUE matlab stub file. not safe to rewrite and use the same path for everything. likely where omni manova run glitches came from
        mat_script = "{}/run_from_python_{}.m".format(self.bash_wrapper_omni_manova_dir, timestamp)
        Path(mat_script).touch()
        Path(log_file).touch()
        completion_code = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        with open(mat_script, 'w') as f:
            with open(mat_script_template, 'r') as old_f:
                old = old_f.read()
            f.write("close all;\n")
            f.write("clear all;\n")
            f.write("save_location='{}';\n".format(out_dir))
            f.write("data_frame_path='{}';\n".format(data_frame_path))
            f.write("test_criteria={};\n".format(test_criteria))
            f.write("log_file='{}';\n".format(log_file))
            f.write("completion_code='{}';\n".format(completion_code))
            f.write("run_R_analysis='{}';\n".format(run_R_analysis))
            f.write(old)
            f.write("exit;")

        cmd = "\"run('{}'); exit;\"".format(mat_script)
        cmd = "{} -nosplash -nodisplay -nodesktop -r {} -logfile {}".format(self.matlab, cmd, log_file)
        if CLUSTER:
            cmd = self.make_cluster_command(cmd, self.bash_wrapper_omni_manova_dir, "omni_manova-{}".format(gen), "120G")
            print(cmd)
            proc = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE)
            jid = proc.stdout.read().decode("utf-8").split(" ")[-1]
            jid = jid.strip("\n")
            cluster_wait_cmd = "pipeline_utilities cluster_wait {}".format(jid)
            print("waiting for Omni Manova call to complete")
            subprocess.run(cluster_wait_cmd.split(" "))
        else:
            subprocess.run(cmd.split(" "))

    
    def make_dataframe(self, generation):
        output_dir = "{}/Omni_Manova-{}".format(self.omni_manova_dir_base, generation)
        #date = datetime.today().strftime('%Y-%m-%d')
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        csv_table = pd.read_csv(self.experiment_table_path, header=0, delimiter=",")
        experiment_list = csv_table.to_dict(orient='records')

        # must use a different template on the cluster because there are many absolute paths in the dataframe
        # unsure how/if possible to use environment variables in that context
        df_template_path = "{}/other/{}_dataframe_template.csv".format(os.path.dirname(os.path.realpath(__file__)), self.project_code)
        if CLUSTER:
             df_template_path = "{}/other/{}_dataframe_template_cluster.csv".format(os.path.dirname(os.path.realpath(__file__)), self.project_code)
        reader = csv.DictReader(open(df_template_path))
        df_template = {}
        for row in reader:
            runno = row["specimen"]
            df_template[runno] = row

        print(type(df_template))
        print(df_template)
        print(df_template["N59128NLSAM"])

        result_csv_dict = {}
        index = 0
        runno_list = df_template.keys()
        debug = self.debug
        quit_index = 0
        # we only want to take experiments from the current (nth) and the immediately previous (n-1) generation to run omni-manova on
        # and we know (by our definition) that each generation has 24 experiments
        if generation >0:
          generation_start_index = SOL_PER_GENERATION * generation
        else:
          generation_start_index = 0

        for experiment in experiment_list:
            if experiment["uid"] < generation_start_index:
                continue
            for runno in runno_list:
                # row is now a dictionary
                # force a COPY here so you do not edit the original
                # and to force each row to be distincy instead of a reference to the same dict over and over
                row = dict(df_template[runno])

                # YAY, supposedely, we can now use real names again. Let's see if that is true...

                # this works, BUT it preserves the useful column names. can't have that
                # want to rename to groupN or subgroupN
                # in this experiment, GROUP is defined to be all of the specimen invariants (strain,sex,age...)
                # and subgroup is defined as all of the dsi studio tractography parameters
                row.update({x: experiment[x] for x in experiment.keys() - self.exclusion_list})
                print(row)
                # subgroup_index = 1
                # for key in sorted(experiment.keys()-exclusion_list):
                #    row["subgroup{}".format(subgroup_index)] = experiment[key]
                #    #print("% subgroup{} = {}".format(subgroup_index, key))
                #    subgroup_index += 1
                # find the connectome file and add it to dict as "file"
                # for now just hard code to count pass conn
                # TODO: this does not allow it to exist in the MDS space
                experiment_number = experiment["uid"]
                row["group_position"] = experiment_number
                # change these to get different connectome types for comparisons
                # i think that omni manova currently only accepts a single connectome file
                c_value = "count"
                c_type = "pass"
                found_connectivity_file = glob.glob("{}/{}/{}/*_labels.count.pass.connectivity.mat".format(self.output_dir_base, experiment_number, runno))
                if len(found_connectivity_file) == 0:
                    print(
                        "no connectivity found, i think this is the end of currently acquired data experiment:{} runno:{}".format(
                            experiment_number, runno))
                    break
                connectome_file = found_connectivity_file[0].replace("\\", "/")
                print("exp={} - runno={} - connect file found: {}".format(experiment_number, runno, connectome_file))
                row["file"] = connectome_file

                # NOOO this is all wrong
                # insidious error here
                # I keep passing by reference
                # so above, when I add and update fields to row, i am also editing that row within df_template
                # and each entry here is set as a reference to row, which in turn is a reference to df_template[runno]
                # so in the end I have a dict with the correct number of rows, but the same parameters for all
                # rows of the same runno
                result_csv_dict[index] = row
                index += 1
            if debug:
                quit_index += 1
            if quit_index > 3:
                break
        # now I have a template dictionary
        # in output csv, I need one row for each experiment for each specimen
        # each of these rows need to start with the template values
        # then we add group fields for the dsi studio parameters
        # then we add the relevant connectome.mat filepath as "file"
        # TODO: more descriptive name for dataframe?
        output_file = "{}/dataframe.csv".format(output_dir)
        df = pd.DataFrame(result_csv_dict)
        df = df.T
        df.to_csv(output_file, sep=",", index_label="Var1")
        '''return()
        # bad way to write below
        csv_columns = list(result_csv_dict[0].keys())
        csv_columns.insert(0, "Var1")
        print(csv_columns)
        with open(output_file, "w") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=csv_columns)
            for key, value in result_csv_dict.items():
                # dictionaries are ordered by insertion value
                # this recreates the dict (inserting the outer key/index at the front of the dict)
                # print(value)
                # print(result_csv_dict[key])
                # continue

                row = {"Var1": key}
                row.update(value)
                writer.writerow(row)'''

    '''def on_generation(self, ga_instance):
        current_sol = ga_instance.population
        current_gen = ga_instance.generations_completed
        #modify the initial_population csv file
        initial_pop_df = pd.read_csv(self.experiment_table_path, header=0, delimiter=",")
        current_sol_df = initial_pop_df.iloc[:24].copy()  # make a new solution df
        # maybe change to 24, bc uid 23 is doubled
        current_sol_df['uid'] = (24*current_gen) + current_sol_df['uid']
        keys = current_sol_df.columns.tolist()
        solution_keys = keys[3:8] # modify accordingly
        current_sol_df[solution_keys] = current_sol
        exp_list = current_sol_df.to_dict(orient='records')
        concat_df = pd.concat([initial_pop_df, current_sol_df], axis=0)
        concat_df.to_csv(self.experiment_table_path,index=False)

        print("starting DSI Studio for generation {}".format(current_gen))
        self.run_dsi_studio(exp_list)
        self.make_dataframe(current_gen)
        print("starting OmniManova for generation {}".format(current_gen))
        self.run_omnimanova(current_gen)
        print("Completed DSI Studion and OmniManova runs for generation {}".format(current_gen))'''

    def on_generation(self, ga_instance):
      #global last_fitness
      print(f"Generation = {ga_instance.generations_completed}")
      print(f"Fitness    = {ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]}")
      best_sol_idx = ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[2]
      self.best_sol_uids.append(ga_instance.generations_completed*SOL_PER_GENERATION + best_sol_idx)
      #print(f"Change     = {ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1] - last_fitness}")
      #last_fitness = ga_instance.best_solution(pop_fitness=ga_instance.last_generation_fitness)[1]



    def on_mutation(self, ga_instance, offspring_mutation):
        """this function should take in the mutated_children
        add them to the run_dsi)studio parameter csv
        run dsi studio
        update omni manova dataframe
        run omni manova
        !! currently is only a direct copy of on_generation, except with mutated_children argument added"""
        current_sol = offspring_mutation
        current_gen = ga_instance.generations_completed + 1
        # modify the initial_population csv file
        initial_pop_df = pd.read_csv(self.experiment_table_path, header=0, delimiter=",")
        # copies the previous generation into a new dataframe that will represent the newly created solution space (previous gen plus current gen to be created)
        current_sol_df = initial_pop_df.iloc[:SOL_PER_GENERATION].copy()  # make a new solution df
        # this increments up the UID of the previous generation, not what we need to do.
        # the previous generation should keep their numbers, and the new generation should be thusly named
        current_sol_df['uid'] = (SOL_PER_GENERATION*(current_gen)) + current_sol_df['uid']
        keys = current_sol_df.columns.tolist()
        solution_keys = keys[3:8]  # modify accordingly
        current_sol_df[solution_keys] = current_sol
        exp_list = current_sol_df.to_dict(orient='records')
        concat_df = pd.concat([initial_pop_df, current_sol_df], axis=0)
        concat_df.to_csv(self.experiment_table_path, index=False)

        print("starting DSI Studio for generation {}".format(current_gen))
        self.run_dsi_studio(exp_list)
        print("making dataframe for generation {}".format(current_gen))
        self.make_dataframe(current_gen)
        print("finished dataframe for generation {}".format(current_gen))
        print("starting OmniManova for generation {}".format(current_gen))
        self.run_omnimanova(current_gen)
        print("Completed DSI Studio and OmniManova runs for generation {}".format(current_gen))

    def on_start(self, ga_instance):
        # code will search for improvements to this by checking what is already completed
        ga_instance.generations_completed = 0
        # this should check where we currently are and try to restart the algorithm from that point
        # do a "generation search" by seeing how many Omni-Manovas have been done
        found_omni_results = glob.glob("{}/Omni_Manova-*".format(self.omni_manova_dir_base))
        # important to sort because the result is not guaranteed to come back already sorted
        found_omni_results.sort()
        omni_done = False
        omni_gens_completed = 0
        last_omni_gen_completed = 0
        for omni_result in reversed(found_omni_results):
            omni_result = omni_result.replace("\\", "/")
            print("searching for omni_results in: {}".format(omni_result))
            # check to make sure that it finished by looking for the Semipar file. this is needed to calculate fitness
            found_semipar = glob.glob("{}/BrainScaled_Omni_Manova/*/Global_Semipar_Image_0000.mat".format(omni_result))
            if len(found_semipar) > 0:
                omni_done = True
                # calculate what generation this is from the omni result path
                last_omni_gen_completed = omni_result.split("-")[-1]
                omni_gens_completed = int(last_omni_gen_completed) + 1
                print("last omni manova gen completed  = {}".format(last_omni_gen_completed))
                print("omni generations completed (count): {}".format(omni_gens_completed))
                break
        if not omni_done:
            # TODO: setup so that it will just run omni manova for the experiments that already have dsi studio results
            # this is a BAD state, as you should already have one generation of manual omni manova results to prime
            # the algorithm. is required for the first round of fitness calculations,
            print("no omni manova results found. BAD STATE. you need a manual first generation from which we calculate initial fitness values for")
            print("search_dir was {}".format(self.omni_manova_dir_base))
            self.make_dataframe(0)
            self.run_omnimanova(0)
            #quit()

        # we are checking if everything that is in our current experiment list CSV file is accounted for.
        # do we need to run DSI studio again before running on_mutation?
        print("checking for dsi_studio results completion")
        # get the list of experiments to check for (similar to on_mutation())
        initial_pop_df = pd.read_csv(self.experiment_table_path, header=0, delimiter=",")
        experiment_list = initial_pop_df.to_dict(orient='records')

        # check that the results are completed (similar to run_dsi_studio())
        dsi_studio_complete = True
        count = 0
        for experiment in experiment_list:
            experiment_dir = "{}/{}".format(self.output_dir_base, experiment["uid"])
            if not os.path.exists(experiment_dir):
                # then the results are not complete yet
                print("DSI Studio results not complete for experiment {}".format(experiment["uid"]))
                dsi_studio_complete = False
                break
            for runno in self.runno_list:
                output_dir = "{}/{}".format(experiment_dir, runno)
                if not os.path.exists(output_dir):
                    # then the results are not complete yet
                    print("DSI Studio results not complete for experiment {} runno {}".format(experiment["uid"], runno))
                    dsi_studio_complete = False
                    break
                # check for completion by looking for any trk file in the outputs folder
                found_trk = glob.glob("{}/*.trk.gz".format(output_dir))
                found_tt = glob.glob("{}/*.tt.gz".format(output_dir))
                if len(found_trk) > 0 or len(found_tt) > 0:
                    # work is completed. congrats
                    count += 1
                else:
                    # then the results are not complete yet
                    print("DSI Studio results not complete for experiment {} runno {}".format(experiment["uid"], runno))
                    dsi_studio_complete = False
                    break

        # try to run dsi studio again, then check for completion after that
        if not dsi_studio_complete:
            # then try to run it
            self.run_dsi_studio(experiment_list)
            # check for completion a second time
            dsi_studio_complete = True
            count = 0
            for experiment in experiment_list:
                experiment_dir = "{}/{}".format(self.output_dir_base, experiment["uid"])
                if not os.path.exists(experiment_dir):
                    # then the results are not complete yet
                    print("DSI Studio results not complete for experiment {}".format(experiment["uid"]))
                    dsi_studio_complete = False
                    break
                for runno in self.runno_list:
                    output_dir = "{}/{}".format(experiment_dir, runno)
                    if not os.path.exists(output_dir):
                        # then the results are not complete yet
                        print("DSI Studio results not complete for experiment {} runno {}".format(experiment["uid"], runno))
                        dsi_studio_complete = False
                        break
                    # check for completion by looking for any trk file in the outputs folder
                    found_trk = glob.glob("{}/*.trk.gz".format(output_dir))
                    found_tt = glob.glob("{}/*.tt.gz".format(output_dir))
                    if len(found_trk) > 0 or len(found_tt) > 0:
                        # work is completed. congrats
                        count += 1
                    else:
                        # then the results are not complete yet
                        print("DSI Studio results not complete for experiment {} runno {}".format(experiment["uid"], runno))
                        dsi_studio_complete = False
                        break


        
        # count how many generations were completed by DSI Studio/
        dsi_gens_completed = 0
        if dsi_studio_complete:
            print("all DSI studio results complete for the existing inputs in {}".format(self.experiment_table_path))
            # count how many generations were completed by DSI Studio/the experiment list
            dsi_gens_completed = len(experiment_list) / SOL_PER_GENERATION
            print("DSI generations completed (count): {}".format(dsi_gens_completed))
        else:
            print("DSI Studio results not complete, but some have been done")
            dsi_gens_completed = count / SOL_PER_GENERATION / len(self.runno_list)
            print("completed {} out of expected {} full generations".format(dsi_gens_completed, len(experiment_list)/SOL_PER_GENERATION))
            # convert to int to truncate it, we don't care about the fractionally complete gen, must redo (will only regenerate incompleted ones)
        dsi_gens_completed = int(dsi_gens_completed)


        # now with the given information about what is complete, make decision about where to start
        if omni_gens_completed < dsi_gens_completed:
            print("DSI results complete, but omni results are beind one generation. running omni manova-{} to catch up, then will continue the main algorithm".format(dsi_gens_completed))
            # TODO: this currently makes the assumption that this case will always be:
                # omni_gens_completed = dsi_gens_completed - 1
            # then we need to catch up by running make dataframe and omni manova
            # we can then update the generation and start the algorithm
            self.make_dataframe(dsi_gens_completed - 1)
            self.run_omnimanova(dsi_gens_completed - 1)
            # TODO: make better names because the indexing here is kind of confusing
            # confusion was that the initial generation 0 doesn't count (bc was manually made and not made by the genetic algo)
            ga_instance.generations_completed = dsi_gens_completed - 1

        if omni_gens_completed == dsi_gens_completed:
            # then no pre-work needs to be done, just set the current generation and go on to main algorithm
            print("parity with omni manova and DSI results. setting generations_completed to {} and resuming main algorithm".format(dsi_gens_completed-1))
            ga_instance.generations_completed = dsi_gens_completed - 1

    def fitness_function(self, ga_instance, solution, solution_idx):
        gen = ga_instance.generations_completed
        # read in global MDS results as a pandas dataframe
        # it will have cryptic column names
        # want to redefine it
        Semipar_path = "{}/Omni_Manova-{}/BrainScaled_Omni_Manova/*/Global_Semipar_Image_0000.mat".format(self.omni_manova_dir_base, gen)
        print(Semipar_path)
        Semipar_path = glob.glob(Semipar_path)
        Semipar_path = Semipar_path[0]
        print("FOUND PATH: {}".format(Semipar_path))
        file = scio.loadmat(Semipar_path) # read semipar file
        cleaned_semipar = [str(item[0]).replace('[', '').replace(']', '').replace("'", "").split() for item in
                        file['full_group_name']]
        semipar_df = pd.DataFrame(cleaned_semipar)
        print(semipar_df)
        # 17 is length of runno list
        # TODO: make a variable for this. this will break future projects
        uids = [(i // 17) for i in range(len(cleaned_semipar))]
        semipar_df['uid'] = uids
        print(semipar_df)
        # TODO: very uncertrain about what gen should be set to, so subtracting one because I get index outofbounds
        if gen == 0:
            # then we only have SOL_PER_GENERATION in total
            dist_data = semipar_df[semipar_df['uid'] == solution_idx]
        else:
            # then we have 2*SOL_PER_GENERATION, but we only care about the final n
            dist_data = semipar_df[semipar_df['uid'] == SOL_PER_GENERATION + solution_idx]
        print("DIST DATA")
        print(dist_data)
        print("uid indices: {}".format(SOL_PER_GENERATION+solution_idx) )
        print("sols: {}".format(SOL_PER_GENERATION))
        print("sol idx {}".format(solution_idx))
        nTg_idx = dist_data[dist_data[0] == 'nTg'].index
        Tg_idx = dist_data[dist_data[0] == 'Tg'].index
        summ = 0
        for i in range(len(Tg_idx)):
            for j in range(len(nTg_idx)):
                dist = file['Dist'][i][j]
                summ += dist
        print("tg_tdx: {}".format(Tg_idx))
        print("ntg_tdx: {}".format(nTg_idx))
        fitness1 = summ / (len(Tg_idx) * len(nTg_idx))
        return fitness1
        #return [fitness1, fitness2]


    def run_GA(self):
        in_pop_df = pd.read_csv(self.experiment_table_path, header=0, delimiter=",").iloc[-SOL_PER_GENERATION:]
        ini_pop = in_pop_df[['fa_threshold', 'turning_angle', 'step_size', 'min_length', 'max_length']].to_numpy()
        ga_instance = pygad.GA(num_generations=self.ngen,
                               num_parents_mating=self.num_parents_mating,
                               num_genes=self.num_genes,
                               fitness_func=self.fitness_function,
                               on_start = self.on_start,
                               parent_selection_type="sus",
                               crossover_type="scattered",
                               gene_space=self.gene_space,
                               gene_type=self.gene_type,
                               on_generation=self.on_generation,
                               on_mutation=self.on_mutation,
                               initial_population=ini_pop,
                               mutation_percent_genes=20,
                               mutation_probability = 0.3,
                               mutation_type="random",
                               keep_elitism = 0,
                               save_best_solutions = True,
                               keep_parents=0)
        ga_instance.run()
        best_solutions = ga_instance.best_solutions
	best_sol_df = pd.DataFrame(best_solutions, columns = ['fa_threshold', 'turning_angle', 'step_size', 'min_length', 'max_length'])
        best_sol_df['Fitness'] = ga_instance.best_solutions_fitness
	best_sol_df['UID'] = self.best_sol_uids
        best_sol_df.to_csv("{}/{}/{}/best_solutions.csv".format(os.environ["BIGGUS_DISKUS"], project_code, project_folder_name))




# experiment organization
# BIGGUS_DISKUS is your working space on any computer
    # cluster = /private_shares/${username}
    # citrix = B:/ProjectSpace/${username}
# within that, make a project_code folder. this includes all data and experiments related to that project. inputs will be here (fib, src, and labels dir)
# within that, make a project_folder folder. This is used to store all results and intermediate data
# within that, are all of the folders that we create. 

if __name__ == "__main__":
    project_code = "20.5xfad.01"
    project_folder_name = "debug_test_10-by-25"
    # number of generations to run?? why did it run 5 instead of 1?
    ngen = 25
    num_parents_mating = 5
    runno_list = ['N59128NLSAM', 'N59130NLSAM', 'N59132NLSAM', 'N59134NLSAM', 'N60076NLSAM', 'N60145NLSAM',
                  'N60149NLSAM', 'N60151NLSAM', 'N60153NLSAM', 'N60155NLSAM', 'N60165NLSAM', 'N60171NLSAM',
                  'N60206NLSAM', 'N60208NLSAM', 'N60213NLSAM', 'N60215NLSAM', 'N60727NLSAM']
    debug = False
    # citrix computer dependent
    if CLUSTER: 
        dsi_studio = "/cm/shared/workstation/aux/dsi_studio_2024-08-14/dsi_studio" 
    else:
        dsi_studio = "//pwp-civm-ctx01/K/CIVM_APPS/dsi_studio_64/dsi_studio_win_cpu_v2024-08-14/dsi_studio.exe"

    experiment_table_path = "{}/{}/{}/genetic_initial_population.csv".format(os.environ["BIGGUS_DISKUS"], project_code, project_folder_name"{}/{}/{}/genetic_initial_population.csv".format(os.environ["BIGGUS_DISKUS"], project_code, project_folder_name))

    # these are the experiment_table columns that we do not care about
    exclusion_list = ["random_seed", "export", "connectivity_threshold", "connectivity_type", "connectivity_value",
                      "threshold_index", "thread_count", "interpolation", "initial_dir", "source", "output",
                      "connectivity", "connectivity_output"]

    # fa_threshold [0,0.7]
    # min_length [0.1, 5]
    # max_len [10,200]
    # turning_angle [15, 60]
    # step_size [0.01, 0.05] step size is in mm. our data has a resolution of 0.025 mm. ranges from sub-voxel to 2-voxels
    # otsu_threshold [0.3, 1.2]
    gene_space = [{'low':1, 'high':70}, [15, 25, 35, 45], {'low':0.01, 'high':0.05},
                  {'low':0.1, 'high':5}, {'low':10, 'high':200}]

    gene_type = [int, int, float, float, int]

    new_search = GA_pipeline(ngen=ngen, experiment_table_path=experiment_table_path, runno_list=runno_list, debug=False,
                             dsi_studio=dsi_studio, exclusion_list=exclusion_list, num_parents_mating=num_parents_mating,
                             gene_space=gene_space, gene_type=gene_type, project_code=project_code)


    new_search.run_GA()

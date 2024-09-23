import pygad
import subprocess
import pandas as pd
import os
import numpy as np
import nibabel as nib
import glob
import json


# value boundaries
# fa_threshold [0,0.7]
# min_length [0.1, 5]
# max_len [10,200]
# turning_angle [15, 60]
# step_size [0.01, 0.05] step size is in mm. our data has a resolution of 0.025 mm. ranges from sub-voxel to 2-voxels
# otsu_threshold [0.3, 1.2]


class GA_pipeline:
    def __init__(self, ngen, project_code, experiment_table_path, output_dir_base,
                 runno_list, debug, dsi_studio, exclusion_list, num_parents_mating, sol_per_pop, gene_space):

        self.project_code = project_code
        self.debug = debug
        self.dsi_studio_version = dsi_studio
        self.exclusion_list = exclusion_list
        self.runno_list = runno_list
        self.experiment_table_path = experiment_table_path
        self.output_dir_base = "B:/ProjectSpace/vc144/{}/{}".format(self.project_code, "parameter_sets")
        self.archive_dir = "A:/{}/research".format(project_code)
        self.src_dir = "B:/ProjectSpace/vc144/{}/{}".format(self.project_code, "src")
        self.fib_dir = "B:/ProjectSpace/vc144/{}/{}".format(self.project_code, "fib")
        self.ngen = ngen
        self.num_parents_mating = num_parents_mating
        self.sol_per_pop = sol_per_pop
        self.num_genes = 12
        self.gene_space = gene_space
        self.matlab = "C:/CIVM_Apps/MATLAB/R2021b/bin/matlab.exe"




    def setup_dsi_studio_trk_call(self, experiment: dict, fib_file, qa_nifti_file, label_file, output_dir):
        options_str = ""
        for key in experiment:
            if key in self.exclusion_list:
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

        s = "{} --action=trk {}".format(DSI_STUDIO, options_str)
        # print(s)
        return s

    def find_vol_pct_threshold(self, image, percent):
        nifti = nib.load(image)
        data, header = nifti.get_fdata(), nifti.header
        data = data[data != 0]
        data = np.sort(data)
        index = round((percent / 100) * len(data))
        threshold = data[index]
        print(threshold)
        return threshold

    def run_dsi_studio(self, sol):
        csv_table = pd.read_csv(self.experiment_table_path, header=1, delimiter="\t")
        keys = csv_table.columns.tolist()
        vals = csv_table.iloc[0].tolist()
        const_keys = list(set(keys) - set(keys[1:13]))
        solution_keys = keys[1:13]
        const_vals = vals[0].append(vals[13:])

        experiment_list = dict(zip(const_keys, const_vals))
        sol_dic = dict(zip(solution_keys, sol))
        experiment_list.update(sol_dic)
        count = 0
        for experiment in [experiment_list]:
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
                label_file = "{}/connectome{}dsi_studio/labels/RCCF/{}_RCCF_labels.nii.gz".format(self.archive_dir, runno,
                                                                                                  runno)
                fib_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz".format(self.fib_dir, runno)
                qa_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz.qa.nii.gz".format(self.fib_dir, runno)

                # check for completion by looking for any trk file in the outputs folder
                found_trk = glob.glob("{}/*.trk.gz".format(output_dir))
                found_tt = glob.glob("{}/*.tt.gz".format(output_dir))
                if len(found_trk) > 0 or len(found_tt) > 0:
                    print("work already complete for experiment {} runno {}".format(experiment["uid"], runno))
                    continue

                cmd = self.setup_dsi_studio_trk_call(experiment, fib_file, qa_file, label_file, output_dir)
                subprocess.run(cmd, shell=True)
            with open('{}/experiment.json'.format(experiment_dir), 'w') as fp:
                json.dump(experiment, fp)
            exit_codes = [p.wait() for p in process_list]

    def run_omnimanova(self, out_dir, data_frame_path):

        mat_script_template = "C:/workstation/code/analysis/Omni_Manova/Run_File/template_prototype_run_from_python.m"
        mat_script = "C:/workstation/code/analysis/Omni_Manova/Run_File/prototype_run_from_python.m"
        test_criteria = "{{'group1','group2','subgroup02','subgroup03', 'subgroup04', 'subgroup05','subgroup06','subgroup07','subgroup08','subgroup09','subgroup10','subgroup11','subgroup12'}}"
        with open(mat_script, 'w') as f:
            with open(mat_script_template, 'r') as old_f:
                old = old_f.read()
            f.write("close all;\n")
            f.write("clear all;\n")
            f.write("save_location='{}';\n".format(out_dir))
            f.write("data_frame_path='{}';\n".format(data_frame_path))
            f.write("test_criteria={};\n".format(test_criteria))
            f.write(old)

        log_file = "{}/omni_manova_{}.log".format(out_dir, time.time())
        cmd = "run('{}'); exit;".format(mat_script)

        cmd = "{} -nosplash -nodisplay -nodesktop -r {} -logfile {}".format(self.matlab, cmd, log_file)

        print(cmd)
        p = subprocess.Popen(cmd)
        if os.path.exists("path_to_global_MDS")==False:
        # this waiting does not work
        #p.wait()
        print("process complete")

    def on_generation(self, ga_instance):
        current_sol = ga_instance.population
        current_gen = ga_instance.generations_completed
        print("starting DSI Studio for generation {}".format(current_gen))
        self.run_dsi_studio(current_sol)
        print("starting OmniManova for generation {}".format(current_gen))
        self.run_omnimanova()
        print("Completed DSI Studion and OmniManova runs for generation {}".format(current_gen))


    def fitness_function(self, ga_instance, solution, solution_idx):
        MDS_data = pd.read_csv("path/to/mds_global_or_ASE") # read MDS global data
        MDS_data = MDS_data[MDS_data['subgroup01'] == solution_idx]
        nTg1 = MDS_data[MDS_data['group1'] == 'nTg'][['X1', 'X2']]
        #nTg2 = MDS_data[MDS_data['group1'] == 'nTg2'][['X1','X2']] # third group
        Tg = MDS_data[MDS_data['group1'] == 'Tg'][['X1', 'X2']]

        Tg_mean = (Tg['X1'].mean(), Tg['X2'].mean()) # average out the MDS points
        nTg1_mean = (nTg1['X1'].mean(), nTg1['X2'].mean()
        #nTg2_mean = (nTg2['X1'].mean(), nTg2['X2'].mean())

        fitness1 = math.dist(Tg_mean, nTg1_mean) # maximize the distance between nTg and Tg groups
        #fitness2 = 1 / math.distance(nTg1_mean, nTg2_mean) # minimize the distance between nTg and agilent groups
        return fitness1
        #return [fitness1, fitness2]

    def run_GA(self):

        ga_instance = pygad.GA(num_generations=self.ngen,
                               sol_per_pop=self.sol_per_pop,
                               num_parents_mating=self.num_parents_mating,
                               num_genes=self.num_genes,
                               fitness_func=fitness_function,
                               parent_selection_type="nsga2",
                               crossover_type="scattered",
                               gene_space=self.gene_space,
                               on_generation=self.on_generation,
                               # initial population
                               mutation_type="adaptive")
        ga_instance.run()
        ga_instance.plot_fitness()





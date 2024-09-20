from datetime import datetime
import os
import pandas as pd
import csv
import glob


# I should have a template/master spreadsheet that keeps the invariant information such as spec id, age, gene_status,etc
# then, for each experiment for each specimen I should create a row in the result spreadsgeet
# this is all fields from the master, plus all relevant fields from the params_csv, plus path to the connectome.mat


if __name__ == '__main__':
    project_code = "20.5xfad.01"
    date = datetime.today().strftime('%Y-%m-%d')
    # the columns we DO NOT want to end up in the omni manova dataframe
    exclusion_list = ["uid", "output", "source", "connectivity", "connectivity_output", "export", "connectivity_value"]
    experiment_table_path = "B:/ProjectSpace/vc144/connectome_parameter_search/connectome_parameter_lut.csv"
    output_dir = "B:/ProjectSpace/vc144/{}/Omni_Manova-{}".format(project_code, date)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    csv_table = pd.read_csv(experiment_table_path, header=1, delimiter="\t")
    experiment_list = csv_table.to_dict(orient='records')

    df_template_path = "B:/ProjectSpace/vc144/connectome_parameter_search/other/20.5xfad.01_dataframe_template.csv"

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
    debug = 0
    if debug:
        runno_list = ["N59128NLSAM", "N59130NLSAM"]
    quit_index = 0
    for experiment in experiment_list:
        quit_index += 0
        for runno in runno_list:
            # row is now a dictionary
            # force a COPY here so you do not edit the original
            # and to force each row to be distincy instead of a reference to the same dict over and over
            row = dict(df_template[runno])
            # this works, BUT it preserves the useful column names. can't have that
            # want to rename to groupN or subgroupN
            # in this experiment, GROUP is defined to be all of the specimen invariants (strain,sex,age...)
            # and subgroup is defined as all of the dsi studio tractography parameters
            #row.update({x: experiment[x] for x in experiment.keys()-exclusion_list})
            #print(row)
            subgroup_index = 1
            for key in sorted(experiment.keys()-exclusion_list):
                row["subgroup{}".format(subgroup_index)] = experiment[key]
                #print("% subgroup{} = {}".format(subgroup_index, key))
                subgroup_index += 1
            # find the connectome file and add it to dict as "file"
            # for now just hard code to count pass conn
            # TODO: this does not allow it to exist in the MDS space
            experiment_number = experiment["uid"]
            row["group_position"] = experiment_number
            # change these to get different connectome types for comparisons
            # i think that omni manova currently only accepts a single connectome file
            c_value = "count"
            c_type = "pass"
            found_connectivity_file = glob.glob("B:/ProjectSpace/vc144/{}/parameter_sets/{}/{}/*_labels.count.pass.connectivity.mat".format(project_code, experiment_number, runno))
            if len(found_connectivity_file) == 0:
                print("no connectivity found, i think this is the end of currently acquired data experiment:{} runno:{}".format(experiment_number, runno))
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
        if quit_index > 5:
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
    exit()
    # bad way to write below
    csv_columns = list(result_csv_dict[0].keys())
    csv_columns.insert(0, "Var1")
    print(csv_columns)
    with open(output_file, "w") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=csv_columns)
        for key, value in result_csv_dict.items():
            # dictionaries are ordered by insertion value
            # this recreates the dict (inserting the outer key/index at the front of the dict)
            #print(value)
            #print(result_csv_dict[key])
            #continue

            row = {"Var1": key}
            row.update(value)
            writer.writerow(row)
            #print(row)

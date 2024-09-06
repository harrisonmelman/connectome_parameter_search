"""
I think it will be much easier to create this pipeline in python
we can use dictionaries/etc to more easily keep organized compared to bash
this is a pseudocode/outline of what we want to get done

we will use the newest version of DSI studio, as we have found that this does not give significant differences from
the legacy pipeline version

src and fib creation is already complete. these files will not vary throughout our experiments
"""
import subprocess
import pandas as pd
import os
import numpy as np
import nibabel as nib
DSI_STUDIO = "K:/CIVM_APPS/dsi_studio_64/dsi_studio_win_cpu_v2024-08-14/dsi_studio.exe"


def main():
    project_code = "20.5xfad.01"
    # need this to grab the label set
    archive_dir = "A:/{}/research".format(project_code)

    experiment_table_path = "B:/ProjectSpace/vc144/connectome_parameter_search/connectome_parameter_lut.csv"

    # in this directory, we will create a sub folder for each
    output_dir_base = "B:/ProjectSpace/vc144/{}/{}".format(project_code, "parameter_sets")

    # I have moved the src and fib files to their own directories for easy reuse
    src_dir = "B:/ProjectSpace/vc144/{}/{}".format(project_code, "src")
    fib_dir = "B:/ProjectSpace/vc144/{}/{}".format(project_code, "fib")

    csv_table = pd.read_csv(experiment_table_path, header=1, delimiter="\t")
    experiment_list = csv_table.to_dict(orient='records')

    runno_list = ["N59128NLSAM"]
    for experiment in experiment_list:
        # load in the parameter set from CSV file
        # an experiment is one line of the CSV
        experiment_dir = "{}/{}".format(output_dir_base, experiment["uid"])
        if not os.path.exists(experiment_dir):
            os.makedirs(experiment_dir)
        for runno in runno_list:
            output_dir = "{}/{}".format(experiment_dir, runno)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            label_file = "{}/connectome{}dsi_studio/labels/RCCF/{}_RCCF_labels.nii.gz".format(archive_dir, runno, runno)
            fib_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz".format(fib_dir, runno)
            qa_file = "{}/nii4D_{}.src.gz.gqi.0.9.fib.gz.qa.nii.gz".format(fib_dir, runno)

            # formulate dsi studio call
            setup_dsi_studio_trk_call(experiment, fib_file, qa_file, label_file, output_dir)

            # use subprocess to run dsi studio
            pass
        pass


def find_vol_pct_threshold(image, percent):
    nifti = nib.load(image)
    data, header = nifti.get_fdata(), nifti.header
    data = data[data != 0]
    data = np.sort(data)
    index = round((percent / 100) * len(data))
    threshold = data[index]
    print(threshold)
    return threshold



# TODO: this does not yet handle FILEPATHS
def setup_dsi_studio_trk_call(experiment: dict, fib_file, qa_nifti_file, label_file, output_dir):
    """takes in a set of DSI_Studio action=trk paramters and returns the command string"""
    options_str = ""
    exclusion_list = ["uid", "fa_threshold", "output", "source", "connectivity", "connectivity_output"]
    for key in experiment:
        if key in exclusion_list:
            continue

        options_str += "--{}={} ".format(key, experiment[key])
    # inputs and outputs, these are different for every individual, so handled separately from the options dict
    options_str += "--output={} ".format(output_dir)
    options_str += "--source={} ".format(fib_file)
    options_str += "--connectivity={} ".format(label_file)
    # if we do not set this, then it will be saved right next to the trk file. we want this.
    #options_str += "--connectivity_output={} ".format(output_dir)

    threshold = find_vol_pct_threshold(qa_nifti_file, experiment["fa_threshold"])
    options_str += "--fa_threshold={} ".format(threshold)

    s = "{} {}".format(DSI_STUDIO, options_str)
    print(s)
    return s


if __name__ == "__main__":
    main()

    # create omni manova "dataframe"
    # call omni manova code to run (MATLAB)

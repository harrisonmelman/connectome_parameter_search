import subprocess
import time
import random
import string
from pathlib import Path
import os

MATLAB = "C:/CIVM_Apps/MATLAB/R2021b/bin/matlab.exe"

out_dir = "B:/ProjectSpace/vc144/20.5xfad.01/omni_manova_python_completion_test"
data_frame_path = "{}/dataframe.csv".format(out_dir)

mat_script_template = "C:/workstation/code/analysis/Omni_Manova/Run_File/template_prototype_run_from_python.m"
mat_script = "C:/workstation/code/analysis/Omni_Manova/Run_File/prototype_run_from_python.m"

# now update the start script with the correct outpath and dataframe path
    # save_location='B:\ProjectSpace\vc144\20.5xfad.01\Omni_Manova-2024-09-19';
    # data_frame_path='B:\ProjectSpace\vc144\20.5xfad.01\Omni_Manova-2024-09-19\dataframe.csv';
# we will also need to update the "test_criteria" variable
    # test_criteria={{'group1','group2','subgroup3', 'subgroup8', 'subgroup18'}}; %5 Sources of variation
# must write to a new file. the only thing we can do to an existing file is to append to the end of it
# which parameters have variation in them?
test_criteria = "{{'group1','group2','subgroup03', 'subgroup04', 'subgroup09'}}"
timestamp = time.time()
log_file = "{}/omni_manova_{}.log".format(out_dir, timestamp)
Path(log_file).touch()
N = 10
completion_code = ''.join(random.choices(string.ascii_uppercase + string.digits, k=N))
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
    f.write(old)


cmd = "run('{}'); exit;".format(mat_script)

cmd = "{} -nosplash -nodisplay -nodesktop -r {} -logfile {}".format(MATLAB, cmd, log_file)

print(cmd)
p = subprocess.Popen(cmd)

# wait until we read completion_code in the log file
with open(log_file, 'r') as f:
    while True:
        time.sleep(5)
        if completion_code in f.read():
            print("FOUND THE CODE: {}".format(completion_code))
            break
print("process complete")

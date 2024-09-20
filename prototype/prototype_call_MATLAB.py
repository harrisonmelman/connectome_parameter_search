# this looks at an approach using subprocess to open git bash with our full environment so that we can call matlab
# no much luck

import subprocess
from subprocess import Popen, PIPE
import nibabel as nib
import numpy as np
#matlab_run -D100 --purpose=tmp_hdr_correct_N58211NLSAM_RCCF_labels --dir_work="/k/ProjectSpace/jjc29/Projects/22.gaj.49/samba_packages/22.gaj.49_pack_tforms" "nii_header_transfer('/k/ProjectSpace/jjc29/Projects/22.gaj.49/22.gaj.49/research/diffusionN58211NLSAMdsi_studio/nii4D_N58211NLSAM.nii','/k/ProjectSpace/jjc29/Projects/22.gaj.49/samba_packages/22.gaj.49_pack_tforms/tmp_hdr_correct_N58211NLSAM_RCCF_labels.nii');"
matlab_run = "K:/workstation/code/shared/pipeline_utilities/matlab_run.pl"
matlab_function = "find_vol_pct_threshold"
pct = 10
img = "B:/ProjectSpace/vc144/20.5xfad.01/fib/nii4D_N59128NLSAM.src.gz.gqi.0.9.fib.gz.qa.nii.gz"
dir_work = "B:/ProjectSpace/vc144/test-pct"
#purpose = "find {}% threshold for {}".format(pct, img)

#cmd = "{} -d100 --purpose='{}' --dir_work='{}' '{}('{}','{}')'".format(matlab_run, purpose, dir_work, matlab_function, img, pct)
#print(cmd)

#cmd = ["/usr/bin/bash", "perl", matlab_run, "-d100", "--purpose={}".format(purpose), "--dir_work={}".format(dir_work), "find_vol_pct_threshold({},{})".format(img,pct)]
#subprocess.run(cmd, shell=True)


process = Popen(["K:/DevApps/Git/git-bash.exe"], stdin=PIPE, stdout=PIPE)
#process.stdin.write(b'Hello\n')
#process.stdin.flush()
# repr attempts to give a print-representable version of the object
#print(repr(process.stdout.readline())) # Should print 'Hello\n'
#process.stdin.write(b'World\n')
#process.stdin.flush()
#print(repr(process.stdout.readline())) # Should print 'World\n'
#process.stdin.close()
print('Waiting for cat to exit')
process.wait()
print('cat finished with return code %d' % process.returncode)


exit()
# calls that I tried that did not work
from subprocess import Popen, PIPE
p = Popen("""   K:/DevApps/Git/usr/bin/echo.exe HOWDY; K:/DevApps/Git/usr/bin/sleep.exe 5;""", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("""K:/DevApps/Git/usr/bin/echo.exe HOWDY; K:/DevApps/Git/usr/bin/sleep.exe 5;""", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("""K:/DevApps/Git/git-bash.exe; K:/DevApps/Git/usr/bin/echo.exe HOWDY; K:/DevApps/Git/usr/bin/sleep.exe 5;""", stdin=PIPE, stdout=PIPE, shell=True)
p = Popen("""K:/DevApps/Git/git-bash.exe; K:/DevApps/Git/usr/bin/echo.exe HOWDY; K:/DevApps/Git/usr/bin/sleep.exe 5;""", stdin=PIPE, stdout=PIPE, shell=True)
p = Popen("""K:/DevApps/Git/git-bash.exe; K:/DevApps/Git/usr/bin/echo.exe HOWDY; K:/DevApps/Git/usr/bin/sleep.exe 5;""", shell=True)
p = Popen("""K:/DevApps/Git/git-bash.exe K:/DevApps/Git/usr/bin/echo.exe HOWDY""", shell=True)
p = Popen("K:/DevApps/Git/git-bash.exe K:/DevApps/Git/usr/bin/echo HOWDY", shell=True)
p = Popen("K:/DevApps/Git/usr/bin/echo HOWDY", shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("/k/DevApps/Git/usr/bin/echo HOWDY", shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("/k/DevApps/Git/usr/bin/echo HOWDY", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("echo HOWDY", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("ls -l /b/", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")

p = Popen("matlab_run -D100 --purpose=test_call_matlab --dir_work=/b/ProjectSpace/vc144 find_vol_pct_threshold(B:/ProjectSpace/vc144/20.5xfad.01/fib/nii4D_N59128NLSAM.src.gz.gqi.0.9.fib.gz.qa.nii.gz,10)", stdin=PIPE, stdout=PIPE, shell=True, executable="K:/DevApps/Git/git-bash.exe")
p = Popen("K:/DevApps/Git/git-bash.exe matlab_run -D100 --purpose=test_call_matlab --dir_work=/b/ProjectSpace/vc144 find_vol_pct_threshold(B:/ProjectSpace/vc144/20.5xfad.01/fib/nii4D_N59128NLSAM.src.gz.gqi.0.9.fib.gz.qa.nii.gz,10)", stdin=PIPE, stdout=PIPE, shell=True)


# some more functions they have that I am not sure how to work exactly yet
# subprocess.call
# subprocess.communicate

dsi_studio="//pwp-civm-ctx01/k/CIVM_Apps/dsi_studio_64/dsi_studio_win_cpu_v2024-08-14/dsi_studio.exe"
in_dir="A:/20.5xfad.01/research"
op_dir_base="B:/ProjectSpace/vc144/20.5xfad.01"

#specimen_list="N59128NLSAM N59130NLSAM N59132NLSAM N59134NLSAM N60076NLSAM N60145NLSAM N60149NLSAM"
#specimen_list="N60151NLSAM N60153NLSAM N60155NLSAM N60165NLSAM N60171NLSAM"
specimen_list="N60206NLSAM N60208NLSAM N60213NLSAM N60215NLSAM N60727NLSAM"
for x in $specimen_list;
do
	diffusion_dir="$in_dir/diffusion${x}dsi_studio";
	op_dir="$op_dir_base/${x}"
	#if [  ! -d $src_op  ]; then
	mkdir $op_dir;
	#fi

	connectome_dir="$in_dir/connectome${x}dsi_studio";
	nii4d_file="$diffusion_dir/nii4D_${x}.nii";
	btable_file="$diffusion_dir/b_table.txt";
	label_file="$connectome_dir/labels/RCCF/${x}_RCCF_labels.nii.gz"
	threshold_file_10pct="${connectome_dir}/threshold_at_10pct_nqa.txt";
	threshold=$(cat $threshold_file_10pct);

	# SRC
	#$dsi_studio --action=src --source=$nii4d_file --b_table=$btable_file --output="$op_dir/nii4D_${x}";
	src_file="$op_dir/nii4D_${x}.src.gz"

    # make a mask
    recmask=${op_dir}/nii4D_${x}_recmask.nii.gz;
    ImageMath 3 ${recmask} ReplaceVoxelValue ${label_file} 1 1.84467440737096e+19 1
    # use mask with --mask $mask option

	#Fib Creation
	if [ ${x:0:2} == "N5" ]; then
		$dsi_studio --action=rec --source=$src_file --thread_count=16 --method=4   --param0=0.9 --other_output=ad,fa,md,rd,qa,iso,rdi,nrdi --mask=${recmask} --check_btable=0 --align_acpc=0 --cmd="[Step T2][B-table][flip by]"
	else
		$dsi_studio --action=rec --source=$src_file --thread_count=16 --method=4   --param0=0.9 --other_output=ad,fa,md,rd,qa,iso,rdi,nrdi --mask=${recmask} --check_btable=0 --align_acpc=0 --cmd="[Step T2][B-table][flip bx]+[Step T2][B-table][flip by]"
	fi

	#Connectome Creation
	$dsi_studio --action=trk --source=$op_dir/*.fib.gz --turning_angle=45 --random_seed=0 --output=$op_dir/*.tt.gz --connectivity=$label_file --min_length=0.5 --export=tdi,tdi_color,stat --seed_plan=0 --connectivity_value=count,ncount,dti_fa,ad,md,rd,mean_length --threshold_index=qa --connectivity_threshold=0 --fiber_count=2000000 --smoothing=0.01 --max_length=200.0 --method=0 --step_size=0.01 --fa_threshold=$threshold --initial_dir=2 --interpolation=0 --connectivity_type=pass
	echo done
done

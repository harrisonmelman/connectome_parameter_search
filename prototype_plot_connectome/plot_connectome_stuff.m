v2022=load('A:\20.5xfad.01\research\connectomeN59130NLSAMdsi_studio\nii4D_N59130NLSAM.src.gqi.0.9.fib.gz.2000K.trk.gz.N59130NLSAM_RCCF_labels.count.pass.connectivity.mat');
v2024=load('B:\ProjectSpace\vc144\20.5xfad.01\bada_btable\N59130NLSAM\nii4D_N59130NLSAM.src.gz.gqi.0.9.tt.gz.N59130NLSAM_RCCF_labels.count.pass.connectivity.mat');

% how to plot a full connectome
figure; title('full connectome comparison N59130')
subplot(1,2,1); imagesc(v2022.connectivity); title('N59130 v2022');
subplot(1,2,2); imagesc(v2024.connectivity); title('N59130 v2022');

row=152;
figure; title('one most significant connectome row N59130');
d = [v2022.connectivity(row,:); v2024.connectivity(row,:)];
stem(d');
legend('v2022', 'v2024');
title('N59130 row 152 stem plot');
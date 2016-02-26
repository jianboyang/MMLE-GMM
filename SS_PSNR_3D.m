function [snr] = SS_PSNR_3D(original, noisy, max_value)

% Calculates PSNR between original image and noisy version

% Hongcheng Wang @UTRC
term = max(original(:));
original = original/term;
noisy = noisy/term;

max_value = max(max(max(abs(original))));
[height,width,nf]=size(original);
error = abs(original - noisy);
enorm = reshape(sum(sum(error.^2)),1,nf);
snr = 10*log10((max_value*max_value*height*width)./enorm);


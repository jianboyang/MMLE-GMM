function [PSNR, SSIM] = saveResults(Xpre, Row, Col, M, T,y,Xtst,filename,Algorim,time,model,nAlgorithm,C)
Xrecon = zeros(Row,Col,M*T); 
PSNR = SS_PSNR_3D(Xtst,Xpre);
Xtst = Xtst/max(Xtst(:));
for k=1:(M*T)
    Xrecon(:,:,k) = Xpre(:,:,k)/max(Xpre(:));
    SSIM(k) = ssim(Xrecon(:,:,k), Xtst(:,:,k));
    MSE(k) = norm(Xrecon(:,:,k)- Xtst(:,:,k),'fro')/norm(Xtst(:,:,k),'fro');
end

if nAlgorithm>30
    savename = [filename '_T' num2str(T) '_F' num2str(M) '_C' num2str(C) '_' Algorim];
else
    savename = [filename '_T' num2str(T) '_F' num2str(M) '_' Algorim];
end
format short
save(savename, '-v7.3', 'Xrecon','Xpre','PSNR','time','SSIM','model','MSE');

writerObj = VideoWriter([savename '.mp4'],'MPEG-4');
writerObj.FrameRate = 12;
open(writerObj);
scrsz = get(0,'ScreenSize');
fig=figure('Position',[50 100 floor(scrsz(3)*0.8) floor(scrsz(4)*0.6)]);
for nF=1:(M*T)
    subplot(1,3,1);
    imshow(Xtst(:,:,nF)./max(Xtst(:)));
    title(['Original Video, Frame:' num2str(nF)]);
    
    
    subplot(1,3,2);
    imshow(y(:,:,ceil(nF/T))/max(max(y(:,:,ceil(nF/T)))));
    title(['Measurement video, Frame: ' num2str(ceil(nF/T))]);
    
    subplot(1,3,3);
    imshow(Xrecon(:,:,nF)./max(Xrecon(:)));
    title(['Recon video, PSNR: ' num2str(PSNR(nF))]);
        
    %pause(0.02);
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);
PSNR = mean(PSNR);
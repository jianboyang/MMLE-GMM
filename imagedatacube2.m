function imagedatacube2(data,mymax,wavelengths,labelss,textposition,cols,rows,writefile,grayc)

%data = datacube to plot
%mymax = scaling factor to increase brightness of each slice
%wavelengths = wavelengths of each spectral slice in nm
%textposition = where the wavelengths text will be placed

n = size(data,1);
% data=mymax*data;
temp=sort(data(:));
temp=mean(temp(end-50:end));
data=data*(50/temp);
if mymax > 1
    data=data*mymax;
end

totalfigs = size(data,3);

% cols = ceil(sqrt(totalfigs));
% rows = ceil(totalfigs/cols);

figscount = 1;
set(gcf,'color','white');

subplot1(rows,cols,'Gap',[0.005 0.005]);
for r = 1:rows
    for c = 1:cols
        if figscount>totalfigs
            subplot1(figscount);axis off;
        else
            currentwl = wavelengths(figscount);

            cmM=(gray*kron(ones(3,1),spectrumRGB(wavelengths(figscount))));
            cmM=cmM/max(max(cmM));

            if figscount<=totalfigs
                subplot1((r-1)*cols+c);
                if isempty(labelss)
                    subimage2(squeeze(data(:,:,figscount)),colormap(cmM));
                    text(textposition(1),textposition(2),['\bf' num2str(currentwl,4) ' nm'],'color','w')
                else
                    if grayc==1
                        subimage2(make01(squeeze(data(:,:,figscount))));
                    else
                        subimage2(squeeze(data(:,:,figscount)),colormap(cmM));
                    end
                    if r*c==1
                        text(textposition(1),textposition(2),['\bf' num2str(labelss(r*c),4) ' frame'],'color','w')
                    else
                        text(textposition(1),textposition(2),['\bf' num2str(labelss((r-1)*cols+c),4) ' frames'],'color','w')
                    end
                end
                if writefile
                    imwrite(squeeze(data(:,:,figscount)),colormap(cmM),['file' num2str(figscount) '.png'])
                end
                axis off;
            end
        end
        figscount = figscount+1;
    end
end

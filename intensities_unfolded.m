% visualize all the data in unolded space
nbins = 1000;
thickbins = 10;
binsize = 1/thickbins;

for LR = 'lr'
    
    load(sprintf('unfolded.%s/indexed/AP_laplace.mat',LR))
    load(sprintf('unfolded.%s/indexed/indexes_mapping.mat',LR))
    load(sprintf('unfolded.%s/indexed/IO_isovolume.mat',LR))
    load(sprintf('unfolded.%s/indexed/PD_laplace.mat',LR))
    mkdir(sprintf('unfolded.%s/visualization',LR));
    orig = load_untouch_nii(sprintf('hippocampus.%s.40um.nii.gz',LR));
    igm = orig.img(idxgm);
    idg = orig.img(idxdg);
    clear orig %save memory
    
    APind = AP_laplace*nbins; PDind = PD_laplace*nbins;
    flatmapped_intensity = nan([nbins nbins thickbins]);
    
    m=0;
    for n = 1:-binsize:binsize
        m=m+1;
        ikeep = find(IO_isovolume<n & IO_isovolume>n-binsize);
        
%         figure('units','normalized','outerposition',[0 0 1 1]);
%         if LR == 'l'
%             title('left hippocampus');
%         elseif LR == 'r'
%             title('right hippocampus');
%         end
%         scatter3(AP_laplace(ikeep),PD_laplace(ikeep),IO_isovolume(ikeep),12,igm(ikeep),'.');
%         axis([0 1 0 1 0 1]);
%         view(-4,75);
%         xlabel('anterior-posterior'); ylabel('proximal-distal'); zlabel('laminar')
%         colormap(jet); h = colorbar; ylabel(h,'staining intensity'); caxis([min(igm) max(igm)]);
%         % Write to the GIF File
%         drawnow;
%         frame = getframe(gcf);
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         if m == 1
%             imwrite(imind,cm,sprintf('unfolded.%s/visualization/stainingIntensity.gif',LR),'gif', 'Loopcount',inf);
%         else
%             imwrite(imind,cm,sprintf('unfolded.%s/visualization/stainingIntensity.gif',LR),'gif','WriteMode','append');
%         end
%         close
        
        %bin
        for AP = 1:nbins
            for PD = 1:nbins
                vox = find((APind<=AP & APind>AP-1) & (PDind<=PD & PDind>PD-1));
                if ~isempty(vox)
                    flatmapped_intensity(AP,PD,m) = mean(igm(vox));
                end
            end
        end
    save(sprintf('unfolded.%s/flatmaps/intensities.mat',LR),'flatmapped_intensity');
    end
end
clear
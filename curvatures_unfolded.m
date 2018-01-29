%curvature - uses downsampled versions to save memory
threshold = [-0.5 0.5]; %min and max curvature at a given vertex
nbins = 1000;

for LR = 'lr'
    
    load(sprintf('downsampled/unfolded.%s/indexed/indexes_mapping.mat',LR))
    load(sprintf('downsampled/unfolded.%s/indexed/IO_laplace.mat',LR))
    
    Space3d = zeros(ds_sz);
    Space3d(ds_idxgm) = IO_laplace;
    
    medSurf = isosurface(Space3d,0.5);
%     medSurf = reducepatch(medSurf,0.1); % reduce number of faces
    meancurve = patchcurvature(medSurf,1); %this step takes a long time

    % delete bad values
%     bad = find(meancurve>=threshold(2) | meancurve<=threshold(1));
%     meancurve(bad) = []; medSurf.vertices(bad,:) = [];
    
    load(sprintf('downsampled/unfolded.%s/indexed/PD_laplace.mat',LR))
    load(sprintf('downsampled/unfolded.%s/indexed/AP_laplace.mat',LR))
    
    AP = zeros(ds_sz);
    AP(ds_idxgm) = AP_laplace;
    PD = zeros(ds_sz);
    PD(ds_idxgm) = PD_laplace;
    
    y = round(medSurf.vertices(:,1)); %have to switch x and y because matlab does surfaces stupidly
    x = round(medSurf.vertices(:,2));
    z = round(medSurf.vertices(:,3));
    
    for n = 1:length(z)
    APind(n) = AP(x(n),y(n),z(n));
    PDind(n) = PD(x(n),y(n),z(n));
    end
    
    save(sprintf('unfolded.%s/indexed/curvatures.mat',LR),'APind','PDind','meancurve','medSurf','threshold');
    
    % plot
    figure('units','normalized','outerposition',[0 0 1 1]);
    scatter3(APind,PDind,meancurve,[],meancurve,'.');
    if LR == 'l'
        title('left hippocampus');
    elseif LR == 'r'
        title('right hippocampus');
    end
    view(-4,75);
    axis([0 1 0 1 threshold]);
    xlabel('anterior-posterior'); ylabel('proximal-distal'); zlabel('mean curvature')
    colormap(jet); colorbar; caxis(threshold);
    % Write to the png File
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,sprintf('unfolded.%s/visualization/curvature.png',LR),'png');
    
    %bin
    APind = APind*nbins; PDind = PDind*nbins;
    flatmapped_curvature = nan([nbins nbins]);
    for AP = 1:nbins
        for PD = 1:nbins
            vox = find((APind<=AP & APind>AP-1) & (PDind<=PD & PDind>PD-1));
            if ~isempty(vox)
                flatmapped_curvature(AP,PD) = mean(meancurve(vox));
            end
        end
    end
%     flatmapped_curvature = inpaintn(flatmapped_curvature);
    save(sprintf('unfolded.%s/flatmaps/curvature.mat',LR),'flatmapped_curvature','threshold');
    
end
clear
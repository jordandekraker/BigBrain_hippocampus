%thickness - uses downsampled versions to save memory
%start points are at isocenter. first stream inward, then outward, then sum
%them together and remove outliers

voxelsize = 0.08; %mm
threshold = [0.1 5.0]; %min and max in mm
nbins = 1000;

stepsize = 0.1;
maxvert = threshold(2)/voxelsize/stepsize;
for LR = 'lr'
    
    load(sprintf('downsampled/unfolded.%s/indexed/indexes_mapping.mat',LR))
    load(sprintf('downsampled/unfolded.%s/indexed/IO_laplace.mat',LR))
    
    % srlm to outside
    Space3d = ones(ds_sz);
    Space3d(ds_idxgm) = IO_laplace;
    fv = isosurface(Space3d,0.5);
    startpts = round([fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3)]);
    clear fv
    
    [dx,dy,dz]=gradient(Space3d);
    streams1 = stream3(dx,dy,dz,startpts(:,1),startpts(:,2),startpts(:,3),[stepsize,maxvert]);
    clear dx dy dz Space3d
    
    % outside to srlm
    Space3d = ones(ds_sz);
    Space3d(ds_idxgm) = -IO_laplace+1;
    
    [dx,dy,dz]=gradient(Space3d);
    streams2 = stream3(dx,dy,dz,startpts(:,1),startpts(:,2),startpts(:,3),[stepsize,maxvert]);
    clear dx dy dz Space3d
    
    % measure the streams
    %     pts = [];
    for n = 1:length(startpts)
        streamlengths1(n) = stepsize*voxelsize*length(streams1{n});
        streamlengths2(n) = stepsize*voxelsize*length(streams2{n});
        %     pts = [pts;streams1{n}];
    end
    streamlengths = streamlengths1 + streamlengths2;
    
    % delete bad streams
%     bad = find(streamlengths1>=threshold(2)/2|streamlengths1<=threshold(1)/2 | streamlengths2>=maxthick/2|streamlengths2<=minthick/2);
%     streamlengths(bad) = []; startpts(bad,:) = [];
    
    % map back to unfolded space
    load(sprintf('downsampled/unfolded.%s/indexed/PD_laplace.mat',LR))
    load(sprintf('downsampled/unfolded.%s/indexed/AP_laplace.mat',LR))
    
    AP = zeros(ds_sz);
    AP(ds_idxgm) = AP_laplace;
    PD = zeros(ds_sz);
    PD(ds_idxgm) = PD_laplace;
    AP = permute(AP,[2 1 3]);
    PD = permute(PD,[2 1 3]);
    
    for n = 1:length(startpts)
        APind(n) = AP(startpts(n,1),startpts(n,2),startpts(n,3));
        PDind(n) = PD(startpts(n,1),startpts(n,2),startpts(n,3));
    end
    
    save(sprintf('unfolded.%s/indexed/thickness.mat',LR),'APind','PDind','streamlengths','startpts','streams1','streams2');
    
    % plot
    figure('units','normalized','outerposition',[0 0 1 1]);
    scatter3(APind,PDind,streamlengths,[],streamlengths,'.');
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
    imwrite(imind,cm,sprintf('unfolded.%s/visualization/thickness.png',LR),'png');
    
    %bin
    APind = APind*nbins; PDind = PDind*nbins;
    flatmapped_thickness = nan([nbins nbins]);
    for AP = 1:nbins
        for PD = 1:nbins
            vox = find((APind<=AP & APind>AP-1) & (PDind<=PD & PDind>PD-1));
            if ~isempty(vox)
                flatmapped_thickness(AP,PD) = mean(streamlengths(vox));
            end
        end
    end
%     flatmapped_thickness = inpaintn(flatmapped_thickness);
    save(sprintf('unfolded.%s/flatmaps/thickness.mat',LR),'flatmapped_thickness','threshold');

end
clear
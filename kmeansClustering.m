% feature clustering
nbins = 1000;
thickbins = 10;
binsize = 1/(thickbins);
x = [binsize:binsize:1] - binsize/2; %note: does this overestimate the CofG? should be 0:?:1, length 10?

for LR = 'lr'
    
    load(sprintf('unfolded.%s/flatmaps/curvature.mat',LR))
    flatmapped_curvature(flatmapped_curvature<threshold(1)) = threshold(1);
    flatmapped_curvature(flatmapped_curvature>threshold(2)) = threshold(2);
    load(sprintf('unfolded.%s/flatmaps/thickness.mat',LR))
    flatmapped_thickness(flatmapped_thickness<threshold(1)) = threshold(1);
    flatmapped_thickness(flatmapped_thickness>threshold(2)) = threshold(2);
    load(sprintf('unfolded.%s/flatmaps/intensities.mat',LR))
    
    flatmapped_intensity = inpaintn(flatmapped_intensity);
    flatmapped_curvature = inpaintn(flatmapped_curvature);
    flatmapped_thickness = inpaintn(flatmapped_thickness);
    
    n = 1/nbins;
    [APgrad,PDgrad] = meshgrid(n:n:1,n:n:1);
    
    flatmapped_intensity = reshape(flatmapped_intensity,nbins*nbins,thickbins);
    flatmapped_curvature = reshape(flatmapped_curvature,nbins*nbins,1);
    flatmapped_thickness = reshape(flatmapped_thickness,nbins*nbins,1);
    APgrad = reshape(APgrad,nbins*nbins,1);
    PDgrad = reshape(PDgrad,nbins*nbins,1);
    
    features(:,1) = mean(flatmapped_intensity,2);
    for n = 1:nbins*nbins %loop through these because CofG calculation is weird..
    features(n,2) = sum(x.*flatmapped_intensity(n,:))/sum(flatmapped_intensity(n,:));
    end
    features(:,3) = moment(flatmapped_intensity,2,2);
    features(:,4) = moment(flatmapped_intensity,3,2);
    features(:,5) = moment(flatmapped_intensity,4,2);
    
    flatmapped_intensity = abs(diff(flatmapped_intensity,1,2));
    features(:,6) = mean(flatmapped_intensity,2);
    for n = 1:nbins*nbins %loop through these because CofG calculation is weird..
    features(n,7) = sum(x.*flatmapped_intensity(n,:))/sum(flatmapped_intensity(n,:));
    end
    features(:,8) = moment(flatmapped_intensity,2,2);
    features(:,9) = moment(flatmapped_intensity,3,2);
    features(:,10) = moment(flatmapped_intensity,4,2);
    
    features(:,11) = flatmapped_curvature;
    features(:,12) = flatmapped_thickness;
    features(:,13) = APgrad;
    features(:,14) = PDgrad;
    clear flatmapped_intensity flatmapped_curvature flatmapped_thickness APgrad PDgrad
    
    features = zscore(features,0,1);
    
    for nclusters = 1:20
        clusters(:,nclusters) = kmeans(features,nclusters);
        
        figure; axis([1 nbins 1 nbins]);
        if LR == 'l'
            title(sprintf('left hippocampus, %d clusters',nclusters));
        elseif LR == 'r'
            title(sprintf('right hippocampus, %d clusters',nclusters));
        end
        imagesc(reshape(clusters(:,nclusters),nbins,nbins));
        xlabel('anterior-posterior'); ylabel('proximal-distal');
        colormap(jet);
    end
    save(sprintf('unfolded.%s/flatmaps/clustering.mat',LR),'clusters','features');

end
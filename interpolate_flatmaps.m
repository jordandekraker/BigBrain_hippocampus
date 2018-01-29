% interpolate flatmaps
res = 1000;

for LR = 'lr'
    
    %curvature
    load(sprintf('unfolded.%s/indexed/curvatures.mat',LR))
    bad = find(isnan(APind) | isnan(PDind));% | meancurve'<threshold(1) | meancurve'>threshold(2));
    meancurve(meancurve<threshold(1)) = threshold(1);
    meancurve(meancurve>threshold(2)) = threshold(2);
    APind(bad) = []; PDind(bad) = []; meancurve(bad) = [];
    APind = ceil(APind*res); PDind = ceil(PDind*res);
    
    APind = [APind, ones(1,res), ones(1,res)*res];
    PDind = [PDind, 1:res, 1:res];
    newvaluesA = zeros(res,1); newvaluesP = zeros(res,1);
%     for AP = 1:res
%         try
%         newvaluesA(AP) = meancurve(find(APind==AP,1,'first'));
%         newvaluesP(AP) = meancurve(find(APind==AP,1,'last'));
%         end
%     end    
    meancurve = [meancurve;newvaluesA;newvaluesP];
    
    
    F = scatteredInterpolant(APind',PDind',meancurve,'natural','none');

    flatmapped_curvature = F({1:res,1:res})';
    imagesc(flatmapped_curvature);
    caxis(threshold); colormap(jet);
    
end
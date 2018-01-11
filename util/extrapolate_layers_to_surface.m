function [dists_layers_outcropping] = extrapolate_layers_to_surface(distance,layheights)
% This function takes the picked layers and extrapolates them to the
% surface assuming a constant slope equal to the average slope from the 
% last "num" samples.
%
% Laura Kehrl, University of Washington, 1/12/2017

depthcutoff = 70;
num = 25;

dists_layers_outcropping = zeros(size(layheights,1),3);
dists_layers_outcropping(:,:) = nan;
for i=2:size(layheights,1)
    nonnan = find(~isnan(layheights(i,:)));
    if layheights(i,nonnan(1)) < depthcutoff
        dzs = diff(layheights(i,nonnan(1:num)));
        dxs = diff(distance(nonnan(1:num)));
        slope = mean(dzs/dxs);

        % Location where we lose the layer
        dists_layers_outcropping(i,1) = distance(nonnan(1)); 
        % Height where we lose the layer
        dists_layers_outcropping(i,2) = layheights(i,nonnan(1));
        % Inferred location where the layer will outcrop at the surface
        dists_layers_outcropping(i,3) = -1*layheights(i,nonnan(1))/slope+distance(nonnan(1));
    end
end



end

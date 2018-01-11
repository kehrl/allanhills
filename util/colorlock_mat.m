function output = colorlock_mat(input,colormap_str,colorrange);
% (C) Nick Holschuh - Penn State University - 2016 (Nick.Holschuh@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The inputs are:
%
% input - matrix of values that correspond to specific colors in a colormap
% colormap_str - string, the name of the prescribed colormap
% colorrange - the caxis limits to impose when generating the rgb-matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The outputs are:
%
% output - n x m x 3 matrix containing rgb values for the input matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nan_ind = find(isnan(input));
inf_ind = find(input == Inf);
input(inf_ind) = max(max(input(find(input ~= Inf))));
n_inf_ind = find(input == -Inf);
input(n_inf_ind) = min(min(input(find(input ~= -Inf))));

if exist('colormap_str') == 0
    cmap = colormap;
    if length(get(gca,'Children')) > 0;
        colorrange = caxis;
    end
else
    cmap = eval(colormap_str);
end


if exist('colorrange') == 0
    colorrange = [min(min(input)) max(max(input))];
end
if colorrange == 0
    
end


color_opts = linspace(colorrange(1),colorrange(2),length(cmap(:,1))-1);
dc = color_opts(2)-color_opts(1);
mincol = color_opts(1);
input = input-mincol;

input = round_to(input,dc);
input = round(input/dc);

output = ind2rgb(input,cmap);
t1 = output(:,:,1);
t2 = output(:,:,2);
t3 = output(:,:,3);
t1(nan_ind) = 1;
t2(nan_ind) = 1;
t3(nan_ind) = 1;
output(:,:,1) = t1;
output(:,:,2) = t2;
output(:,:,3) = t3;

end





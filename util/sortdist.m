function [x_out,y_out] = sortdist(x_in,y_in)
% Sort input coordinates x,y by closest distance, starting with the first
% index.
%
% Laura Kehrl, University of Washington, 1/12/2017

N = length(x_in);
x_out = zeros(N,1);
y_out = zeros(N,1);

% Set first coordinate
[~,minind] = min(y_in);
x_out(1) = x_in(minind);
y_out(1) = y_in(minind);

x_in(minind) = NaN;
y_in(minind) = NaN;


for i=2:N
    [dist,minind] = min(sqrt((x_out(i-1)-x_in).^2+(y_out(i-1)-y_in).^2));
    x_out(i) = x_in(minind);
    y_out(i) = y_in(minind);
    
    x_in(minind) = NaN;
    y_in(minind) = NaN;
end

end


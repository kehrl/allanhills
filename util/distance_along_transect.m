function d = distance_along_transect(x,y)
% Finds distance between points along a transect
%
% Laura Kehrl, University of Washington, 12/12/2015


N = length(x);
d = zeros(N,1);
for i=2:N
    d(i) = d(i-1)+sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
end


end


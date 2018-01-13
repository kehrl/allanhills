function [xvec,yvec] = find_north(x,y)
%Find the direction that points north based on polar stereographic
%coordinates
%
% Laura Kehrl, University of Washington, 1/12/2017

[sp(1),sp(2)] = ll2ps(-90,0);
sp_slope = (sp(2)-y)/(sp(1)-x);

xdelta = 1;
ydelta = xdelta*sp_slope;

xvec = xdelta/sqrt(xdelta^2+ydelta^2);
yvec = ydelta/sqrt(ydelta^2+ydelta^2);

end


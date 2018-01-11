function [u,v] = velazimuth2coords(x,y,velmag,azimuth)
%This function takes velocities, given as a magnitude velmag with an azimuth 
%in WGS84, to u,v components. Note that the azimuth is given in deg E.
%We have to do this to deal with Nicky's velocity measurements (ugh).
%
% Laura Kehrl, University of Washington, 1/12/2017

u = zeros(length(x),1);
v = zeros(length(y),1);

for i=1:length(x)
    [xvec_north,yvec_north] = find_north(x(i),y(i));
    alpha_north = atan2(yvec_north,xvec_north);
    alpha_total = (pi/180)*(((azimuth(i))))+pi;
    u(i) = velmag(i)*cos(alpha_total);
    v(i) = velmag(i)*sin(alpha_total);
end


end


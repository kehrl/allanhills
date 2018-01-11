function [h_main, h_inset]=inset(main_handle, inset_handle,xsize,ysize)

% The function plotting figure inside figure (main and inset) from 2 existing figures.
% inset_size is the fraction of inset-figure size, default value is 0.35
% The outputs are the axes-handles of both.
% 
% An examle can found in the file: inset_example.m
% 
% Moshe Lindner, August 2010 (C).

if nargin==2
    xsize=0.35;
    ysize=0.35;
end

new_fig=figure;
units = get(main_handle,'units');
set(gcf,'units',units);
set(gcf,'pos',get(main_handle,'pos'));
main_fig = findobj(main_handle,'Type','axes');
h_main = copyobj(main_fig,new_fig);
inset_fig = findobj(inset_handle,'Type','axes');
h_inset = copyobj(inset_fig,new_fig);
set(0,'currentfigure',new_fig);
ax = get(gca,'pos');
xsize=(ax(3)-ax(1))*xsize;
ysize=(ax(4)-ax(2))*ysize;
set(h_inset,'pos', [ax(1)+ax(3)-xsize-0.03*(ax(3)-ax(1)) ax(2)+ax(4)-ysize-0.03*(ax(4)-ax(2)) xsize ysize])
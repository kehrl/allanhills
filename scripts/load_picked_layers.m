% An example script for loading a picked layer from spicker.
%
% Laura Kehrl, UW, 2/13/2017

layers_track1 = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_1b_layers_max.splayer'));

figure(1)
imagesc(track1B.filtdata1B_all,[-0.0015,0.0015]); hold on;
colormap('bone')
for i = 1:size(layers_track1,1)

    plot(1:size(layers_track1,2),layers_track1(i,:),'r','linewidth',1.5);
    
end
set(gca,'ydir','reverse')
xlim([0,size(layers_track1,2)])
ylim([0,size(track1B.filtdata1B_all,1)])
title('Track 1')

layers_track4 = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata4_layers_max.splayer'));

figure(2)
imagesc(track4.filtdata4_all,[-0.0015,0.0015]); hold on;
colormap('bone')
for i = 1:size(layers_track4,1)

    plot(layers_track4(i,:),'r','linewidth',1.5);
    
end
set(gca,'ydir','reverse')
xlim([0,size(layers_track4,2)])
ylim([0,size(track4.filtdata4_all,1)])
title('Track 4')

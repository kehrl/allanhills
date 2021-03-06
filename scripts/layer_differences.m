% This script plots layers in the radar data to look at changes in layer
% thickness with distance. 
%
% Laura Kehrl, UW, 10/01/2017


% Options
deviation_from_mean = 1;


%% Track 1B

% Get desired layers from track1B
n = 1;
for i=1:size(track1B.layages,1)
    nonnan = find(~(isnan(track1B.layheights(i,:))));
    if (track1B.Distance_total1B_all(nonnan(1)) < 5.2e3) && (track1B.Distance_total1B_all(nonnan(end)) > 18e3)
        layers(n,:) = track1B.Height1B_all - track1B.layheights(i,:)';
        layer_ages(n,:) = track1B.layages(i);
        layer_ages_errors(n,:) = track1B.layerrors(i);
        n=n+1;
    end
end
layers = layers(1:end-1,:);
layer_ages = layer_ages(1:end-1);
[~,ind] = sort(layers(:,2000));
layers = flipud(layers(ind,:));
layer_ages = flipud(layer_ages(ind));

track1B.layerpacks = layers;
track1B.layerpacks_ages = layer_ages;


% Calculate layer differences
for i=1:size(layers,1)-1
    layerdiffs(i,:) = layers(i,:)-layers(i+1,:);
    if deviation_from_mean   
        layerdiffs(i,:) = layerdiffs(i,:)-nanmean(layerdiffs(i,:));
    else
        layerdiffs(i,:) = layerdiffs(i,:);
    end
end
track1B.layerpacks_diffs = layerdiffs;

zdata = 800:1:2200;
distdata = track1B.Distance_total1B_all;
filtdata = zeros(length(zdata),length(track1B.Distance_total1B_all));
filtdata(:,:) = nan;
ind = find(track1B.z_out ~= 0);
for i=1:length(track1B.Distance_total1B_all)
    [~,j] = min(abs(track1B.Height1B_all(i)-track1B.thick(i)-zdata));
    filtdata(1:end,i) = interp1(track1B.Height1B_all(i)-track1B.z_out(ind-1:end),track1B.filtdata1B_all(ind-1:end,i),zdata(1:end)); 
end

% Get surface slope and curvature
x_temp = fliplr([track1B.Distance_total1B_all(1):2.5:track1B.Distance_total1B_all(end)]);
zs_temp = interp1(fliplr(track1B.Distance_total1B_all),track1B.Height1B_all,x_temp);

[b,a]=butter(4,0.0011,'low');
zs_temp(1:6) = zs_temp(7);
zs_temp_mean = filtfilt(b,a,zs_temp);

dhdx = gradient(zs_temp,x_temp);
dhdx_mean = filtfilt(b,a,dhdx);
d2hdx2_mean = (gradient(dhdx_mean,x_temp));


%% Plot 1B figure
figure('units','inches');
coloptions = [[0,0,1];[0,1,1];[212/255,175/255,55/255];[1,0,0];[0,0,0]];

hold off;
pos = get(gcf,'pos');
set(gca,'XTickLabel','','YTickLabel','','box','off','visible','off');
set(gcf,'pos',[0 0 3.5 3],'color','white');

[ha, pos] = tight_subplot(2,1,[0.01 0.01],[0.25 0.02],[0.25 0.01]);

subplot(ha(2));
hold on;
box('on');
if deviation_from_mean
    h(1) = plot(fliplr(track1B.Distance_total1B_all/1e3)+dist_offset/1e3,track1B.Height1B_all-layers(1,:)'-nanmean(track1B.Height1B_all-layers(1,:)'),'color',coloptions(1,:),'linewidth',3);
else
    h(1) = plot(fliplr(track1B.Distance_total1B_all/1e3)+dist_offset/1e3,track1B.Height1B_all-layers(1,:)','color',coloptions(1,:),'linewidth',3); 
end
for i=1:size(layerdiffs,1)
    plot(fliplr(track1B.Distance_total1B_all/1e3+dist_offset/1e3),layerdiffs(i,:)','color',coloptions(i+1,:),'linewidth',3);
end
xlim([dist_offset/1e3,track1B.Distance_total1B_all(end)/1e3+dist_offset/1e3])
%ylim([-100,100])
if deviation_from_mean
    ylabel({'Layer thickness';'anomaly (m)'},'fontsize',8,'fontname','Arial')
    ylim([-120,120])
    set(ha(2),'ytick',[-100:50:100],'fontsize',8,'fontname','Arial')
else
    ylabel('Layer thickness (m)','fontsize',8,'fontname','Arial')
    ylim([0,300])
    set(ha(2),'ytick',[0:50:300],'fontsize',8,'fontname','Arial')
end
set(gca,'ticklength',[0.02 0.05])
text(0.5+dist_offset/1e3,100,'c','fontsize',9,'fontname','arial','fontweight','bold');
xlabel('{\it{x}} (km)','fontsize',8,'fontname','Arial');
set(gca,'fontsize',8,'fontname','Arial');
clear h

subplot(ha(1));
pcolor(fliplr(distdata)/1e3+dist_offset/1e3,zdata,filtdata);
shading interp;
set(gca,'clim',[-0.0015 0.0015]);
colormap(bone);
%xlabel('Distance (km)','fontsize',10,'fontname','Arial')
ylabel({'Elevation (m)'},'fontsize',8,'fontname','Arial');
freezeColors;
hold on;
nonnan = find(~(isnan(layers(1,:))));
x_area1 = fliplr(track1B.Distance_total1B_all)+dist_offset/1e3;
x_area = [x_area1(nonnan)/1e3+dist_offset/1e3,fliplr(x_area1(nonnan)/1e3)+dist_offset/1e3];
y_area = [(track1B.Height1B_all(nonnan))',fliplr(layers(1,nonnan))];
h(1) = fill(x_area,y_area,coloptions(1,:),'facealpha',0.5); 
for i=1:size(layers,1)-1
    nonnan = find((~(isnan(layers(i,:)))) & (~(isnan(layers(i+1,:)))));
    x_area = [x_area1(nonnan)/1e3+dist_offset/1e3,fliplr(x_area1(nonnan)/1e3)+dist_offset/1e3];
    y_area = [(layers(i,nonnan)),fliplr(layers(i+1,nonnan))];
    h(i+1) = fill(x_area,y_area,coloptions(i+1,:),'facealpha',0.5,'edgecolor','none');  
end
plot(fliplr(track1B.Distance_total1B_all/1e3)+dist_offset/1e3,track1B.Height1B_all,'k','linewidth',1.25); hold on;
plot(fliplr(track1B.Distance_total1B_all/1e3)+dist_offset/1e3,layers,'k','linewidth',1.25);
plot(fliplr(track1B.Distance_total1B_all/1e3)+dist_offset/1e3,(track1B.Height1B_all-track1B.thick),'k','linewidth',1.5);
set(gca,'ydir','normal')
ylim([800,2200])
xlim([dist_offset/1e3,track1B.Distance_total1B_all(end)/1e3+dist_offset/1e3])
set(ha(1),'XTickLabel','');
set(ha(1),'layer','top','fontsize',8,'fontname','Arial');
set(gca,'ticklength',[0.02 0.05])
text(0.5+dist_offset/1e3,1440,'1''','fontsize',8,'fontname','arial','fontweight','bold');
text(18.9+dist_offset/1e3,1440,'1 ','fontsize',8,'fontname','arial','fontweight','bold');
text(0.5+dist_offset/1e3,940,'a','fontweight','bold','fontsize',9,'fontname','arial');
% Legend
rectangle('Position',[13.05+dist_offset/1e3,875,6.35,455],'facecolor','w','edgecolor','k');
rectangle('Position',[13.3+dist_offset/1e3,1185,0.6,80],'facecolor',coloptions(1,:),'edgecolor','k');
text(14.2+dist_offset/1e3,1240,'1','fontsize',8,'fontname','arial');
rectangle('Position',[14.9+dist_offset/1e3,1185,0.6,80],'facecolor',coloptions(4,:),'edgecolor','k');
text(15.8+dist_offset/1e3,1240,sprintf('4 (%d ka)',round(track1B.layerpacks_ages(4))),'fontsize',8,'fontname','arial');
rectangle('Position',[13.3+dist_offset/1e3,1065,0.6,80],'facecolor',coloptions(2,:),'edgecolor','k');
text(14.2+dist_offset/1e3,1110,'2','fontsize',8,'fontname','arial');
rectangle('Position',[14.9+dist_offset/1e3,1065,0.6,80],'facecolor',coloptions(5,:),'edgecolor','k');
text(15.8+dist_offset/1e3,1110,sprintf('5 (%d ka)',round(track1B.layerpacks_ages(5))),'fontsize',8,'fontname','arial');
rectangle('Position',[13.3+dist_offset/1e3,925,0.6,80],'facecolor',coloptions(3,:),'edgecolor','k');
text(14.2+dist_offset/1e3,980,'3','fontsize',8,'fontname','arial');
set(gca,'XTickLabel','');
clear h 

export_fig(fullfile(REPO_HOME,'figures/ahills_track1B_thickness.pdf'),'-opengl','-r600');
close;


%% Track 4 
clear layers layerdiffs

% Get desired layers from track4
n = 1;
for i=1:size(track4.layages,1)
    nonnan = find(~(isnan(track4.layheights(i,:))));
    if (track4.Distance_total4_all(nonnan(1)) < 6e3) && (track4.Distance_total4_all(nonnan(end)) > 13e3)
        layers(n,:) = track4.Height4_all - track4.layheights(i,:)';
        n=n+1;
    end
end
[~,ind]=sort(layers(:,3000));
layers = flipud(layers(ind,:));
layers = layers([1,7,13,20,end],:);

% Calculate layer differences
for i=1:size(layers,1)-1
    layerdiffs(i,:) = layers(i,:)-layers(i+1,:);
    layerdiffs(i,:) = layerdiffs(i,:)-nanmean(layerdiffs(i,:));
end

zdata = 800:1:2200;
distdata = track4.Distance_total4_all;
filtdata = zeros(length(zdata),length(track4.Distance_total4_all));
filtdata(:,:) = nan;
ind = find(track4.z_out ~= 0);
for i=1:length(track4.Distance_total4_all)
    [~,j] = min(abs(track4.Height4_all(i)-track4.thick(i)-zdata));
    filtdata(1:end,i) = interp1(track4.Height4_all(i)-track4.z_out(ind-1:end),track4.filtdata4_all(ind-1:end,i),zdata(1:end)); 
end

% % plot picked layers and bed for track 4
% figure('units','inches');
% 
% hold off;
% pos = get(gcf,'pos');
% set(gca,'XTickLabel','','YTickLabel','','box','off');
% set(gcf,'pos',[0 0 7 4],'color','white')
% [ha, pos] = tight_subplot(2,1,[0.03 0.03],[0.1 0.01],[0.1 0.01]);
% 
% subplot(ha(2));
% hold on;
% box('on');
% h(1) = plot(fliplr(track4.Distance_total4_all/1e3),track4.Height4_all-layers(1,:)'-nanmean(track4.Height4_all-layers(1,:)'),'k','linewidth',1.5);
% h(2:5) = plot(fliplr(track4.Distance_total4_all/1e3),layerdiffs,'linewidth',1.5);
% legend(h, '0-1','1-2','2-3','3-4','4-5','location','southeast')
% xlim([0,track4.Distance_total4_all(end)/1e3])
% ylim([-100,100])
% set(ha(2),'ytick',[-80:40:80],'fontsize',8,'fontname','Arial')
% ylabel('Deviation from mean thickness (m)','fontsize',8,'fontname','Arial')
% xlabel('Distance (km)','fontsize',8,'fontname','Arial')
% 
% subplot(ha(1));
% pcolor(fliplr(distdata)/1e3,zdata,filtdata);
% shading interp;
% set(gca,'clim',[-0.0015 0.0015]);
% colormap(bone);
% %xlabel('Distance (km)','fontsize',10,'fontname','Arial')
% ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
% %freezeColors;
% hold on;
% h(1) = plot(fliplr(track4.Distance_total4_all/1e3),track4.Height4_all,'k','linewidth',1.5);
% h(2:6) = plot(fliplr(track4.Distance_total4_all/1e3),layers,'linewidth',1.5);
% plot(fliplr(track4.Distance_total4_all/1e3),track4.Height4_all-track4.thick,'k:','linewidth',1.5);
% set(gca,'ydir','normal')
% ylim([800,2200])
% xlim([0,track4.Distance_total4_all(end)/1e3])
% set(ha(1),'XTickLabel','');
% set(ha(1),'layer','top','fontsize',8,'fontname','Arial');
% legend_h = columnlegend(2,h,{'0','1','2','3','4','5'},'location','southeast','boxon');
% clear h
% 
% export_fig(fullfile(REPO_HOME,'figures/ahills_track4_thickness.pdf'),'-zbuffer','-p0.01','-nofontswap');
% close;

clear layers layerdiffs


%% Track 3 

% Get desired layers from track3
n = 1;
for i=1:size(track3.layages,1)
    nonnan = find(~(isnan(track3.layheights(i,:))));
    if (track3.Distance_total3_all(nonnan(1)) < 6e3) && (track3.Distance_total3_all(nonnan(end)) > 10e3)
        layers(n,:) = track3.Height3_all - track3.layheights(i,:)';
        n=n+1;
    end
end
[~,ind]=sort(layers(:,3000));
layers = flipud(layers(ind,:));
layers = layers([1,3,4,9,10],:);

% Calculate layer differences
for i=1:size(layers,1)-1
    layerdiffs(i,:) = layers(i,:)-layers(i+1,:);
    layerdiffs(i,:) = layerdiffs(i,:)-nanmean(layerdiffs(i,:));
end

zdata = 800:1:2200;
distdata = track3.Distance_total3_all;
filtdata = zeros(length(zdata),length(track3.Distance_total3_all));
filtdata(:,:) = nan;
ind = find(track3.z_out ~= 0);
for i=1:length(track3.thick)
    [~,j] = min(abs(track3.Height3_all(i)-track3.thick(i)-zdata));
    filtdata(1:end,i) = interp1(track3.Height3_all(i)-track3.z_out(ind-1:end),track3.filtdata3_all(ind-1:end,i),zdata(1:end)); 
end

% % plot picked layers and bed for track 3
% figure('units','inches');
% 
% hold off;
% pos = get(gcf,'pos');
% set(gca,'XTickLabel','','YTickLabel','','box','off');
% set(gcf,'pos',[0 0 7 4],'color','white')
% [ha, pos] = tight_subplot(2,1,[0.03 0.03],[0.1 0.01],[0.1 0.01]);
% 
% subplot(ha(2));
% hold on;
% box('on');
% h(1) = plot(fliplr(track3.Distance_total3_all/1e3),track3.Height3_all-layers(1,:)'-nanmean(track3.Height3_all-layers(1,:)'),'k','linewidth',1.5);
% h(2:5) = plot(fliplr(track3.Distance_total3_all/1e3),layerdiffs,'linewidth',1.5);
% legend(h, '0-1','1-2','2-3','3-4','4-5','location','southeast')
% xlim([0,track3.Distance_total3_all(end)/1e3])
% ylim([-100,100])
% set(ha(2),'ytick',[-80:40:80],'fontsize',8,'fontname','Arial')
% ylabel('Deviation from mean thickness (m)','fontsize',8,'fontname','Arial')
% xlabel('Distance (km)','fontsize',8,'fontname','Arial')
% 
% subplot(ha(1));
% pcolor(fliplr(distdata)/1e3,zdata,filtdata);
% shading interp;
% set(gca,'clim',[-0.0015 0.0015]);
% colormap(bone);
% %xlabel('Distance (km)','fontsize',10,'fontname','Arial')
% ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
% %freezeColors;
% hold on;
% h(1) = plot(fliplr(track3.Distance_total3_all/1e3),track3.Height3_all,'k','linewidth',1.5);
% h(2:6) = plot(fliplr(track3.Distance_total3_all/1e3),layers,'linewidth',1.5);
% plot(fliplr(track3.Distance_total3_all(1:length(track3.thick))/1e3),track3.Height3_all(1:length(track3.thick))-track3.thick,'k:','linewidth',1.5);
% set(gca,'ydir','normal')
% ylim([800,2200])
% xlim([0,track3.Distance_total3_all(end)/1e3])
% set(ha(1),'XTickLabel','');
% set(ha(1),'layer','top','fontsize',8,'fontname','Arial');
% legend_h = columnlegend(2,h,{'0','1','2','3','4','5'},'location','southeast','boxon');
% clear h
% 
% export_fig(fullfile(REPO_HOME,'figures/ahills_track3_thickness.pdf'),'-zbuffer','-p0.01','-nofontswap');
% close;

clear layers layerdiffs


%% Track 2 


% Get desired layers from track2
n = 1;
for i=1:size(track2.layages,1)
    nonnan = find(~(isnan(track2.layheights(i,:))));
    if (track2.Distance_total2_all(nonnan(1)) < 6e3) && (track2.Distance_total2_all(nonnan(end)) > 13e3)
        layers(n,:) = track2.Height2_all - track2.layheights(i,:)';
        n=n+1;
    end
end
[~,ind]=sort(layers(:,3000));
layers = flipud(layers(ind,:));
layers = layers([1,2,3,end-1,end],:);

% Calculate layer differences
for i=1:size(layers,1)-1
    layerdiffs(i,:) = layers(i,:)-layers(i+1,:);
    layerdiffs(i,:) = layerdiffs(i,:)-nanmean(layerdiffs(i,:));
end

zdata = 800:1:2200;
distdata = track2.Distance_total2_all;
filtdata = zeros(length(zdata),length(track2.Distance_total2_all));
filtdata(:,:) = nan;
ind = find(track2.z_out ~= 0);
for i=1:length(track2.thick)
    [~,j] = min(abs(track2.Height2_all(i)-track2.thick(i)-zdata));
    filtdata(1:end,i) = interp1(track2.Height2_all(i)-track2.z_out(ind-1:end),track2.filtdata2_all(ind-1:end,i),zdata(1:end)); 
end

% plot picked layers and bed for track 2
% figure('units','inches');
% 
% hold off;
% pos = get(gcf,'pos');
% set(gca,'XTickLabel','','YTickLabel','','box','off');
% set(gcf,'pos',[0 0 7 4],'color','white')
% [ha, pos] = tight_subplot(2,1,[0.03 0.03],[0.1 0.01],[0.1 0.01]);
% 
% subplot(ha(2));
% hold on;
% box('on');
% h(1) = plot(fliplr(track2.Distance_total2_all/1e3),track2.Height2_all-layers(1,:)'-nanmean(track2.Height2_all-layers(1,:)'),'k','linewidth',1.5);
% h(2:5) = plot(fliplr(track2.Distance_total2_all/1e3),layerdiffs,'linewidth',1.5);
% legend(h, '0-1','1-2','2-3','3-4','4-5','location','southeast')
% xlim([0,track2.Distance_total2_all(end)/1e3])
% ylim([-100,100])
% set(ha(2),'ytick',[-80:40:80],'fontsize',8,'fontname','Arial')
% ylabel('Deviation from mean thickness (m)','fontsize',8,'fontname','Arial')
% xlabel('Distance (km)','fontsize',8,'fontname','Arial')
% 
% subplot(ha(1));
% pcolor(fliplr(distdata)/1e3,zdata,filtdata);
% shading interp;
% set(gca,'clim',[-0.0015 0.0015]);
% colormap(bone);
% %xlabel('Distance (km)','fontsize',10,'fontname','Arial')
% ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
% %freezeColors;
% hold on;
% h(1) = plot(fliplr(track2.Distance_total2_all/1e3),track2.Height2_all,'k','linewidth',1.5);
% h(2:6) = plot(fliplr(track2.Distance_total2_all/1e3),layers,'linewidth',1.5);
% plot(fliplr(track2.Distance_total2_all(1:length(track2.thick))/1e3),track2.Height2_all(1:length(track2.thick))-track2.thick,'k:','linewidth',1.5);
% set(gca,'ydir','normal')
% ylim([800,2200])
% xlim([0,track2.Distance_total2_all(end)/1e3])
% set(ha(1),'XTickLabel','');
% set(ha(1),'layer','top','fontsize',8,'fontname','Arial');
% legend_h = columnlegend(2,h,{'0','1','2','3','4','5'},'location','southeast','boxon');
% clear h
% 
% export_fig(fullfile(REPO_HOME,'figures/ahills_track2_thickness.pdf'),'-zbuffer','-p0.01','-nofontswap');
% close;

clear layers layerdiffs
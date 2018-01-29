% This script considers the 3D internal stratigraphy. 
%
% Laura Kehrl, UW, 01/09/2018

%% Get location of potential core
potential_core_dists = 16e3;

% Get location
[~,ind_core] = min(abs(potential_core_dists-dist));

%% Get depth of layer across radar grid
% Chosen age
[~,ind] = min(abs(track4.layages-105));
age = track4.layages(ind);
cutoff = 1;

% Track 4 
layerdepth = [track4.X4_all,track4.Y4_all,...
    track4.Height4_all,track4.thick,track4.layheights(ind,:)'];

% Track 1B
[diff,ind] = min(abs(track1B.layages-age));
if diff < cutoff
    layerdepth = [layerdepth; track1B.X1B_all,track1B.Y1B_all,...
    track1B.Height1B_all,track1B.thick,track1B.layheights(ind,:)'];
end

% Track 2
[diff,ind] = min(abs(track2.layages-age));
if diff < cutoff
    layerdepth = [layerdepth; track2.X2_all(1:length(track2.thick)),...
        track2.Y2_all(1:length(track2.thick)),...
        track2.Height2_all(1:length(track2.thick)),...
        track2.thick,track2.layheights(ind,(1:length(track2.thick)))'];
end

% Track 3
[diff,ind] = min(abs(track3.layages-age));
if diff < cutoff
    layerdepth = [layerdepth; track3.X3_all(1:length(track3.thick)),...
        track3.Y3_all(1:length(track3.thick)),...
        track3.Height3_all(1:length(track3.thick)),...
        track3.thick,track3.layheights(ind,(1:length(track3.thick)))'];
end

% Track A
[diff,ind] = min(abs(trackA.layages_track4-age));
if diff < cutoff
    layerdepth = [layerdepth; trackA.X(1:length(trackA.thick)),...
        trackA.Y(1:length(trackA.thick)),...
        trackA.Height(1:length(trackA.thick)),...
        trackA.thick,trackA.layheights(ind,(1:length(trackA.thick)))'];
end
    
% Track B
[diff,ind] = min(abs(trackB.layages_track4-age));
if diff < cutoff
    layerdepth = [layerdepth; trackB.X(1:length(trackB.thick)),...
        trackB.Y(1:length(trackB.thick)),...
        trackB.Height(1:length(trackB.thick)),...
        trackB.thick,trackB.layheights(ind,(1:length(trackB.thick)))'];
end

% Track C
[diff,ind] = min(abs(trackC.layages_track4-age));
if diff < cutoff
    layerdepth = [layerdepth; trackC.X(1:length(trackC.thick)),...
        trackC.Y(1:length(trackC.thick)),...
        trackC.Height(1:length(trackC.thick)),...
        trackC.thick,trackC.layheights(ind,(1:length(trackC.thick)))'];
end

% Track D
[diff,ind] = min(abs(trackD.layages_track4-age));
if diff < cutoff
    layerdepth = [layerdepth; trackD.X(1:length(trackD.thick)),...
        trackD.Y(1:length(trackD.thick)),...
        trackD.Height(1:length(trackD.thick)),...
        trackD.thick,trackD.layheights(ind,(1:length(trackD.thick)))'];
end

% Track E
[diff,ind] = min(abs(trackE.layages_track4-age));
if diff < cutoff
    layerdepth = [layerdepth; trackE.X(1:length(trackE.thick)),...
        trackE.Y(1:length(trackE.thick)),...
        trackE.Height(1:length(trackE.thick)),...
        trackE.thick,trackE.layheights(ind,(1:length(trackE.thick)))'];
end

%% Scatter plot

fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 2 3.6],'color','white');

[x,y,Z] = bedmap2_data('bedw',[500000 520000 520000 500000],[-1300000 -1400000 -1400000 -1300000],'xy');
[c,ch] = contour(x/1e3,y/1e3,Z,[1000:100:2000],'k'); hold on;
set(gca,'fontsize',6,'fontname','arial');

scatter(layerdepth(:,1)/1e3,layerdepth(:,2)/1e3,5,layerdepth(:,4)-layerdepth(:,5)); hold on;
set(gca,'clim',[350,650]);
colormap(flipud(cubehelix([],0.6,-0.9,1.9,1.1,[0,0.8],[0,0.9])))
h = colorbar('southoutside');
xlabel(h,{'Height of 105-ka ice','above bed (m)'},'fontsize',8,'fontname','Arial');
set(h,'ticklength',[0.12 0.12],'fontsize',8);
pos = get(h,'position');
set(h,'position',[pos(1)+0.17 pos(2)-0.02 pos(3)-0.19 pos(4)*0.5]);

axis equal
plot(xtrack(ind_core)/1e3,ytrack(ind_core)/1e3,'kp','markerfacecolor','r','markersize',10)
set(gca,'ticklength',[0.02 0.02]);

text(512.5,-1343.8,'d','fontsize',10,'fontweight','bold','fontname','arial');
text(510.6,-1354,'Northing (km)','rotation',90,'fontsize',8,'fontname','arial')
%text(511.6,-1361,'-1360','rotation',90,'fontsize',6,'fontname','arial')
text(511.6,-1356,'-1355','rotation',90,'fontsize',6,'fontname','arial')
text(511.6,-1351,'-1350','rotation',90,'fontsize',6,'fontname','arial')
text(511.6,-1346,'-1345','rotation',90,'fontsize',6,'fontname','arial')

xlim([min(layerdepth(:,1))/1e3,max(layerdepth(:,1))/1e3+1]);
ylim([min(layerdepth(:,2))/1e3+2.5,max(layerdepth(:,2))/1e3-1.5]);
xlabel('Easting (km)','fontsize',8,'fontname','arial');
tl = clabel(c,ch,'manual','fontsize',7,'fontname','arial','labelspacing',500);

pos = get(gca,'Position');
set(gca,'Position',[pos(1)+0.05 pos(2)+0.18 pos(3)+0.15 pos(4)-0.18]);
set(gca,'xtick',[512:2:518])
set(gca,'ytick',[-1360:5:-1340])
set(gca,'yticklabel',[])

export_fig(fullfile(REPO_HOME,'figures/layer_above_bed.pdf'),'-painters','-transparent')
close;

%% Set up radar data for fence plot

zdata = 700:1:2200;
distdata = track1B.Distance_total1B_all;
filtdata1B = zeros(length(zdata),length(track1B.Distance_total1B_all));
filtdata1B(:,:) = nan;
ind = find(track1B.z_out ~= 0);
for i=1:length(track1B.Distance_total1B_all)
    [~,j] = min(abs(track1B.Height1B_all(i)-track1B.thick(i)-zdata));
	filtdata1B(1:end,i) = interp1(track1B.Height1B_all(i)-track1B.z_out(ind-1:end),track1B.filtdata1B_all(ind-1:end,i),zdata(1:end)); 
end

distdata = trackA.Distance;
filtdataA = zeros(length(zdata),length(trackA.Distance));
filtdataA(:,:) = nan;
ind = find(trackA.z_out ~= 0);
for i=1:length(trackA.Distance)
    [~,j] = min(abs(trackA.Height(i)-trackA.thick(i)-zdata));
	filtdataA(1:end,i) = interp1(trackA.Height(i)-trackA.z_out(ind-1:end),trackA.filtdata(ind-1:end,i),zdata(1:end)); 
end

distdata = trackB.Distance;
filtdataB = zeros(length(zdata),length(trackB.Distance));
filtdataB(:,:) = nan;
ind = find(trackB.z_out ~= 0);
for i=1:length(trackB.Distance)
    [~,j] = min(abs(trackB.Height(i)-trackB.thick(i)-zdata));
	filtdataB(1:end,i) = interp1(trackB.Height(i)-trackB.z_out(ind-1:end),trackB.filtdata(ind-1:end,i),zdata(1:end)); 
end

distdata = trackD.Distance;
filtdataD = zeros(length(zdata),length(trackD.Distance));
filtdataD(:,:) = nan;
ind = find(trackD.z_out ~= 0);
for i=1:length(trackD.Distance)
    [~,j] = min(abs(trackD.Height(i)-trackD.thick(i)-zdata));
	filtdataD(1:end,i) = interp1(trackD.Height(i)-trackD.z_out(ind-1:end),trackD.filtdata(ind-1:end,i),zdata(1:end)); 
end

distdata = trackE.Distance;
filtdataE = zeros(length(zdata),length(trackE.Distance));
filtdataE(:,:) = nan;
ind = find(trackE.z_out ~= 0);
for i=1:length(trackE.Distance)
    [~,j] = min(abs(trackE.Height(i)-trackE.thick(i)-zdata));
	filtdataE(1:end,i) = interp1(trackE.Height(i)-trackE.z_out(ind-1:end),trackE.filtdata(ind-1:end,i),zdata(1:end)); 
end

%% Fence plot

fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 5 3.6],'color','white');

fence(trackA.X/1e3,trackA.Y/1e3,zdata,filtdataA,0.2,5);
hold on;
view([-1 -2 0.8])
fence(trackD.X/1e3,trackD.Y/1e3,zdata,filtdataD,0.2,5)
fence(track1B.X1B_all(1500:end-0)/1e3,track1B.Y1B_all(1500:end-0)/1e3,zdata,filtdata1B(:,1500:end-0),0.2,5)
plot3(trackA.X/1e3,trackA.Y/1e3,trackA.Height-trackA.thick,'k','linewidth',3);
plot3(trackD.X/1e3,trackD.Y/1e3,trackD.Height-trackD.thick,'k','linewidth',3);
plot3(track1B.X1B_all(1500:end-0)/1e3,track1B.Y1B_all(1500:end-0)/1e3,track1B.Height1B_all(1500:end-0)-track1B.thick(1500:end-0),'k','linewidth',3);

ages = [89, 105, 187];
colors = [[178 34 249]; [0 0 255]; [0 128 0]]/255; %from colorbrewer2.org
for i=1:length(ages)
    [~,ind] = min(abs(trackA.layages_track4-ages(i)));
    plot3(trackA.X/1e3,trackA.Y/1e3,trackA.Height-trackA.layheights(ind,:)','color',colors(i,:),'linewidth',3);
    [~,ind] = min(abs(trackD.layages_track4-ages(i)));
    plot3(trackD.X/1e3,trackD.Y/1e3,trackD.Height-trackD.layheights(ind,:)','color',colors(i,:),'linewidth',3);
    [~,ind] = min(abs(track1B.layages-ages(i)));
    plot3(track1B.X1B_all(1500:end-0)/1e3,track1B.Y1B_all(1500:end-0)/1e3,track1B.Height1B_all(1500:end-0)-track1B.layheights(ind,(1500:end-0))','color',colors(i,:),'linewidth',3);
end

% Plot grid lines
set(gca,'xgrid','on','ygrid','on')

% Plot axes
plot3([512,520],[-1340,-1340],[2200,2200],'k')
plot3([520,520],[-1356,-1340],[2200,2200],'k')
plot3([520,520],[-1356,-1356],[600,2200],'k')

h = get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[1 1 h(3)/1.5])
set(gca,'clim',[-0.0015,0.0015])
colormap('bone')
set(gca,'xtick',[510:5:520],'xlim',[512 520])
set(gca,'ytick',[-1360:5:-1340],'ylim',[-1356,-1340])
set(gca,'ztick',[800:400:2200])
zlim([600 2200])
set(gca,'fontsize',6,'fontname','arial');

zlabel('Elevation (m)','fontsize',8,'fontname','arial');
text(507,-1348,620,'Northing (km)','fontsize',8,'fontname','arial');
text(513,-1359,620,'Easting (km)','fontsize',8,'fontname','arial');

text(515,-1350,620,'A','fontsize',8,'fontname','arial','fontweight','bold');
text(514.5,-1347.7,620,'D','fontsize',8,'fontname','arial','fontweight','bold');
text(512.8,-1340.5,680,'1','fontsize',8,'fontname','arial','fontweight','bold');

text(517.5,-1355,1820,'  89 ka','fontsize',8,'fontname','arial')
text(517.5,-1355,1700,'105 ka','fontsize',8,'fontname','arial');
text(517.5,-1355,1580,'187 ka','fontsize',8,'fontname','arial');
text(517.5,-1355,1460,'bed','fontsize',8,'fontname','arial');

export_fig(fullfile(REPO_HOME,'figures/fence_plot.pdf'),'-opengl','-transparent','-r600');
close;
% This script loads age constraints, radar, data, and picked layers for 
% the Allan Hills 2016 field season. It also produces a basemap and 
% plots for the picked layers, calculates age uncertainties for the picked
% layers, and does some preliminary Dansgaard Johnsen to determine a
% preliminary age for ice near the bed.
%
% Laura Kehrl, UW, 3/13/2017


%% Load age constraints 
labels_cores = ['BIT-58';'S27   '];
labels_cores_xy = [510427.0493501444, -1354818.0534980388; 
                   512730.29779785156, -1357364.657766041]; 

% Load S27 ages for plotting on Track 4
S27_vertical_core = importdata(fullfile(REPO_HOME,'ages/spaulding_table5_S27.txt'));
    
% Load surface ages from Table 5 in Spaulding et al., 2013
S27_surface_core = importdata(fullfile(REPO_HOME,'ages/spaulding_table5_surface.txt'));
% Table 5 doesn't have coordinates for the surface ages, so we have to use 
% the distances from S1 to interpolate to actual coordinates. I'm using the
% distance vs. coordinates in the Isotopes.xlsx file for interpolation.
S1_dist = importdata(fullfile(REPO_HOME,'ages/distance_S1_flowline.txt')); % from BIT-58 Surface in Isotopes.xlsx
S27_surface_core_loc = [interp1(S1_dist.data(:,1),S1_dist.data(:,2),S27_surface_core.data(:,2)),...
                    interp1(S1_dist.data(:,1),S1_dist.data(:,3),S27_surface_core.data(:,2))];
clear S1_dist
                
% Load big black coordinates
bigblack = importdata(fullfile(REPO_HOME,'ages/bigblack_polar.txt')); %115 ka
[bigblack(:,1),bigblack(:,2)] = sortdist(bigblack(:,1),bigblack(:,2));

% Load BIT 12/16 coordinates
bit16 = importdata('ages/bit16_polar.txt'); %205 ka

% Load stake measurements (velocities, mass balance, emergence velocities,
% etc.) from Spaulding 2012 paper
[stakes,stakes_txt,stakes_raw] = xlsread(fullfile(REPO_HOME,'velocities/velocities_table1_spaulding.xlsx'));
[stakes_u,stakes_v] = velazimuth2coords(stakes(:,1),stakes(:,2),stakes(:,4),stakes(:,6)); 


%% Load radar data

% Across tracks A,B,C,D,E
% All across tracks are in "crosstracks.mat", so we load that file and
% create new structures for each across track for ease of use.

crosstracks = load('crosstracks.mat');
% Across track A
trackA.filtdata = crosstracks.filtdata_a;
trackA.z_out = crosstracks.z_out;
trackA.X = crosstracks.X1a;
trackA.Y = crosstracks.Y1a;
trackA.Distance = distance_along_transect(crosstracks.X1a,crosstracks.Y1a);
trackA.Height = crosstracks.Height1a;

% ACross track B
trackB.filtdata = crosstracks.filtdata_b;
trackB.z_out = crosstracks.z_out;
trackB.X = crosstracks.X1b;
trackB.Y = crosstracks.Y1b;
trackB.Distance = flipud(distance_along_transect(crosstracks.X1b,crosstracks.Y1b));
trackB.Height = crosstracks.Height1b;

% Across track C
trackC.filtdata = crosstracks.filtdata_c;
trackC.z_out = crosstracks.z_out;
trackC.X = crosstracks.X1c;
trackC.Y = crosstracks.Y1c;
trackC.Distance = distance_along_transect(crosstracks.X1c,crosstracks.Y1c);
trackC.Height = crosstracks.Height1c;

% Across track D
trackD.filtdata = crosstracks.filtdata_d;
trackD.z_out = crosstracks.z_out;
trackD.X = crosstracks.X1d;
trackD.Y = crosstracks.Y1d;
trackD.Distance = flipud(distance_along_transect(crosstracks.X1d,crosstracks.Y1d));
trackD.Height = crosstracks.Height1d;

% Across track E
trackE.filtdata = crosstracks.filtdata_e;
trackE.z_out = crosstracks.z_out;
trackE.X = crosstracks.X1e;
trackE.Y = crosstracks.Y1e;
trackE.Distance = distance_along_transect(crosstracks.X1e,crosstracks.Y1e);
trackE.Height = crosstracks.Height1e;

% Across track F
trackF = load('track2to1.mat');

% Tracks along flow
% Each track has its own file, so we just need to load the mat files.

% Track 1
track1B = load('track1B_all.mat');

% Track 2
track2 = load('track2_all.mat');

% Track 3
track3 = load('track3_all.mat');

% Track 4
track4 = load('track4_all.mat');

% Track BIT58_1
BIT58_1 = load('BIT58_1.mat');

% Track BIT58_2
BIT58_2 = load('BIT58_2.mat');


%% Correct tracks for offset between transmitter and receiver
offset = 18; % distance between gps and center of transmitter/receiver

% Track 1B
newdist = track1B.Distance_total1B_all-offset;
track1B.X1B_all = interp1(track1B.Distance_total1B_all,track1B.X1B_all,newdist)';
track1B.Y1B_all = interp1(track1B.Distance_total1B_all,track1B.Y1B_all,newdist)';
track1B.Distance_total1B_all = newdist;

% Track 2
newdist = track2.Distance_total2_all-offset;
track2.X2_all = interp1(track2.Distance_total2_all,track2.X2_all,newdist)';
track2.Y2_all = interp1(track2.Distance_total2_all,track2.Y2_all,newdist)';
track2.Distance_total2_all = newdist;

% Track 3
newdist = track3.Distance_total3_all-offset;
track3.X3_all = interp1(track3.Distance_total3_all,track3.X3_all,newdist)';
track3.Y3_all = interp1(track3.Distance_total3_all,track3.Y3_all,newdist)';
track3.Distance_total3_all = newdist;

% Track 4
newdist = track4.Distance_total4_all-offset;
track4.X4_all = interp1(track4.Distance_total4_all,track4.X4_all,newdist)';
track4.Y4_all = interp1(track4.Distance_total4_all,track4.Y4_all,newdist)';
track4.Distance_total4_all = newdist;

% Track A
newdist = trackA.Distance-offset;
trackA.X = interp1(trackA.Distance,trackA.X,newdist);
trackA.Y = interp1(trackA.Distance,trackA.Y,newdist);
trackA.Distance = newdist;

% Track B
newdist = trackB.Distance-offset;
trackB.X = interp1(trackB.Distance,trackB.X,newdist);
trackB.Y = interp1(trackB.Distance,trackB.Y,newdist);
trackB.Distance = newdist;

% Track C
newdist = trackC.Distance-offset;
trackC.X = interp1(trackC.Distance,trackC.X,newdist);
trackC.Y = interp1(trackC.Distance,trackC.Y,newdist);
trackC.Distance = newdist;

% Track D
newdist = trackD.Distance-offset;
trackD.X = interp1(trackD.Distance,trackD.X,newdist);
trackD.Y = interp1(trackD.Distance,trackD.Y,newdist);
trackD.Distance = newdist;

% Track E
newdist = trackE.Distance-offset;
trackE.X = interp1(trackE.Distance,trackE.X,newdist);
trackE.Y = interp1(trackE.Distance,trackE.Y,newdist);
trackE.Distance = newdist;

% Track F
newdist = trackF.Distance_total-offset;
trackF.X = interp1(trackF.Distance_total,trackF.X,newdist)';
trackF.Y = interp1(trackF.Distance_total,trackF.Y,newdist)';
trackF.Distance_total = newdist';

% Track BIT58_1
newdist = BIT58_1.Distance_total-offset;
BIT58_1.X = interp1(BIT58_1.Distance_total,BIT58_1.X,newdist)';
BIT58_1.Y = interp1(BIT58_1.Distance_total,BIT58_1.Y,newdist)';
BIT58_1.Distance_total = newdist';

% Track BIT58_2
newdist = BIT58_2.Distance_total-offset;
BIT58_2.X = interp1(BIT58_2.Distance_total,BIT58_2.X,newdist)';
BIT58_2.Y = interp1(BIT58_2.Distance_total,BIT58_2.Y,newdist)';
BIT58_2.Distance_total = newdist';

clear offset newdist


%% Find intersection points between radar transects

% Track 1B
track1B.interx_trackA = interx([trackA.X,trackA.Y]',[track1B.X1B_all,track1B.Y1B_all]');
track1B.interx_trackB = interx([trackB.X,trackB.Y]',[track1B.X1B_all,track1B.Y1B_all]');
track1B.interx_trackC = interx([trackC.X,trackC.Y]',[track1B.X1B_all,track1B.Y1B_all]');
track1B.interx_trackD = interx([trackD.X,trackD.Y]',[track1B.X1B_all,track1B.Y1B_all]');
track1B.interx_trackE = interx([trackE.X,trackE.Y]',[track1B.X1B_all,track1B.Y1B_all]');
track1B.interx_trackF = interx([trackF.X,trackF.Y]',[track1B.X1B_all,track1B.Y1B_all]');

trackA.interx_track1B = track1B.interx_trackA;
trackB.interx_track1B = track1B.interx_trackB;
trackC.interx_track1B = track1B.interx_trackC;
trackD.interx_track1B = track1B.interx_trackD;
trackE.interx_track1B = track1B.interx_trackE;
trackF.interx_track1B = track1B.interx_trackF;

% Track 2
track2.interx_trackA = interx([trackA.X,trackA.Y]',[track2.X2_all,track2.Y2_all]');
track2.interx_trackB = interx([trackB.X,trackB.Y]',[track2.X2_all,track2.Y2_all]');
track2.interx_trackC = interx([trackC.X,trackC.Y]',[track2.X2_all,track2.Y2_all]');
track2.interx_trackD = interx([trackD.X,trackD.Y]',[track2.X2_all,track2.Y2_all]');
track2.interx_trackE = interx([trackE.X,trackE.Y]',[track2.X2_all,track2.Y2_all]');
track2.interx_trackF = interx([trackF.X,trackF.Y]',[track2.X2_all,track2.Y2_all]');

trackA.interx_track2 = track2.interx_trackA;
trackB.interx_track2 = track2.interx_trackB;
trackC.interx_track2 = track2.interx_trackC;
trackD.interx_track2 = track2.interx_trackD;
trackE.interx_track2 = track2.interx_trackE;
trackF.interx_track2 = track2.interx_trackF;

%Track 3
track3.interx_trackA = interx([trackA.X,trackA.Y]',[track3.X3_all,track3.Y3_all]');
track3.interx_trackB = interx([trackB.X,trackB.Y]',[track3.X3_all,track3.Y3_all]');
track3.interx_trackC = interx([trackC.X,trackC.Y]',[track3.X3_all,track3.Y3_all]');
track3.interx_trackD = interx([trackD.X,trackD.Y]',[track3.X3_all,track3.Y3_all]');
track3.interx_trackE = interx([trackE.X,trackE.Y]',[track3.X3_all,track3.Y3_all]');

trackA.interx_track3 = track3.interx_trackA;
trackB.interx_track3 = track3.interx_trackB;
trackC.interx_track3 = track3.interx_trackC;
trackD.interx_track3 = track3.interx_trackD;
trackE.interx_track3 = track3.interx_trackE;

% Track 4
track4.interx_trackA = interx([trackA.X,trackA.Y]',[track4.X4_all,track4.Y4_all]');
track4.interx_trackB = interx([trackB.X,trackB.Y]',[track4.X4_all,track4.Y4_all]');
track4.interx_trackC = interx([trackC.X,trackC.Y]',[track4.X4_all,track4.Y4_all]');
track4.interx_trackD = interx([trackD.X,trackD.Y]',[track4.X4_all,track4.Y4_all]');
track4.interx_trackE = interx([trackE.X,trackE.Y]',[track4.X4_all,track4.Y4_all]');

trackA.interx_track4 = track4.interx_trackA;
trackB.interx_track4 = track4.interx_trackB;
trackC.interx_track4 = track4.interx_trackC;
trackD.interx_track4 = track4.interx_trackD;
trackE.interx_track4 = track4.interx_trackE;

  
%% Get age constraints for Track 4
% S27 is along track 4, so this is the transect with the best age
% constraints. There are age constraints from both the horizontal and
% vertical cores. First we load the picked bed and layers, then we find the
% location of the age constraints, then we use the age constraints to
% determine the age of the picked layers.
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata4_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata4_layers_max.splayer'));
thick = interp1(1:length(track4.z_out),track4.z_out,bedpicks)';
bed = [track4.Longitude4_all(1:length(bedpicks)), track4.Latitude4_all(1:length(bedpicks)),... 
       track4.X4_all(1:length(bedpicks)), track4.Y4_all(1:length(bedpicks)),...
       track4.Height4_all(1:length(bedpicks)),thick];
track4.thick = thick;
track4.bed = track4.Height4_all(1:length(bedpicks))-thick;
labels_xy = [track4.X4_all(20)-350, track4.Y4_all(20)-500; 
             track4.X4_all(length(bedpicks))-1600, -1.34468e6];
labels = ['4 ';'4'''];

% Find location of S27 and surface ages along transect 4 
[~,ind_S27_track4] = min(sqrt((track4.X4_all-labels_cores_xy(2,1)).^2+(track4.Y4_all-labels_cores_xy(2,2)).^2));
dist_S27_track4 = track4.Distance_total4_all(ind_S27_track4)/1e3; 
dist_S27_surface_core_track4 = zeros(length(S27_surface_core_loc),1);
for i=1:length(S27_surface_core_loc)
    [mindist,ind_S27_surface_core_track4] = min(sqrt((track4.X4_all-S27_surface_core_loc(i,1)).^2+(track4.Y4_all-S27_surface_core_loc(i,2)).^2));
    disp(sprintf('Mindist between surface age constraints and track 4 is %f m',mindist));
    dist_S27_surface_core_track4(i) = track4.Distance_total4_all(ind_S27_surface_core_track4);
end
clear mindist ind_S27_surface_core_track4

% Find location of big black along transect
bigblack_track4_loc=interx(bigblack',[track4.X4_all,track4.Y4_all]');
[~,ind_bigblack_track4] = min(sqrt((track4.X4_all-bigblack_track4_loc(1)).^2+(track4.Y4_all-bigblack_track4_loc(2)).^2));
dist_bigblack_track4 = track4.Distance_total4_all(ind_bigblack_track4)/1e3;

% Find location of bit16 along transect
bit16_track4_loc=interx(bit16',[track4.X4_all,track4.Y4_all]');
[~,ind_bit16_track4] = min(sqrt((track4.X4_all-bit16_track4_loc(1)).^2+(track4.Y4_all-bit16_track4_loc(2)).^2));
dist_bit16_track4 = track4.Distance_total4_all(ind_bit16_track4)/1e3;

% Find depth of layers
for i=1:size(laypicks,1)
      track4.layheights(i,:) = interp1(track4.z_out,laypicks(i,:))';
end

% We need to extrapolate the picked layers to the surface so we can give
% them an age. Eventually we will use Seth's radar to do this, but in the
% meantime I can just take the mean slope from the last picked part of the layer
% and extrapolate it to the surface.
dists_layers_outcropping_track4 = extrapolate_layers_to_surface(track4.Distance_total4_all,track4.layheights);
dists_layers_outcropping_track4([23,38],:) = nan; %bad slope for last pick, so this layer intersects with other layers before reaching the surface
track4.layages_outcropping_track4 = interp1(dist_S27_surface_core_track4,S27_surface_core.data(:,1),dists_layers_outcropping_track4(:,3));
track4.dists_layers_outcropping_track4 = dists_layers_outcropping_track4;

% Find age of layers that intersect S27
track4.layages_S27 = interp1(S27_vertical_core.data(:,2),S27_vertical_core.data(:,1),track4.layheights(:,ind_S27_track4));
% Some of the layers intersect S27, but are lost near the surface due to
% ringing. In that case, we use the outcropping distance to infer the depth where
% they would intersect S27.
for i = 1:size(laypicks,1)
    if dists_layers_outcropping_track4(i,3) < dist_S27_track4*1e3 && isnan(track4.layages_S27(i))
        z_S27 = interp1([dists_layers_outcropping_track4(i,1),...
            dists_layers_outcropping_track4(i,3)],[dists_layers_outcropping_track4(i,2),0],...
            dist_S27_track4*1e3);
        track4.layages_S27(i) = interp1(S27_vertical_core.data(:,2),S27_vertical_core.data(:,1),z_S27);
    end
end
clear z_S27

track4.layages = nanmean([track4.layages_outcropping_track4,track4.layages_S27],2);

clear ind_S27_surface_core_track4 ind_bigblack_track4 ind_bit16_track4 


%% Plot Track 4
zdata = 700:1:2200;
distdata = track4.Distance_total4_all;
filtdata = zeros(length(zdata),length(track4.Distance_total4_all));
filtdata(:,:) = nan;
ind = find(track4.z_out ~= 0);
for i=1:length(track4.Distance_total4_all)
    [~,j] = min(abs(track4.Height4_all(i)-track4.thick(i)-zdata));
	filtdata(1:end,i) = interp1(track4.Height4_all(i)-track4.z_out(ind-1:end),track4.filtdata4_all(ind-1:end,i),zdata(1:end)); 
end
%ind = find(isnan(filtdata));
%fildata(ind) = -1e3;

% Plot bed and location of age constraints
figure('units','inches');
hold off;
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 4.8 2.3],'color','white')
% Plot radar track4
uimagesc(distdata/1e3,zdata,filtdata,[-0.0015,0.0015]); hold on;
colormap(flipud(colormap('bone')));
set(gca,'ydir','normal');
%uimagesc(track4.Distance_total4_all/1e3,track4.z_out(ind:end),track4.filtdata4_all(ind:end,:),[-0.0015,0.0015]);
colormap(flipud(colormap('bone')));
freezeColors('nancolor',[1 1 1]);
hold on;
% Plot S27
surf = interp1(distdata,track4.Height4_all,dist_S27_track4*1e3);
plot([dist_S27_track4,dist_S27_track4],[surf,surf-max(S27_vertical_core.data(:,2))],'k','linewidth',2);
text(dist_S27_track4-0.4,1600,'S27','fontsize',8,'fontname','Arial','backgroundcolor','w');
%plot([dist_bigblack_track4,dist_bigblack_track4],[-20,60],'b','linewidth',2);
%plot([dist_bit16_track4,dist_bit16_track4],[-20,60],'b','linewidth',2);
% Plot S27 surface ages from Table 5
surf = interp1(distdata,track4.Height4_all,dist_S27_surface_core_track4);
plot(dist_S27_surface_core_track4/1e3,surf,'ko','markerfacecolor','r','markersize',4);
% Legend
rectangle('Position',[0.25,850,3.7,250],'facecolor','w','edgecolor','k');
%plot([0.4,0.7],[780,780],'b','linewidth',2)
%text(0.9,780,'Tephra layer','fontsize',8,'fontname','arial');
plot(0.5,930,'ko','markerfacecolor','r','markersize',4)
text(0.9,930,'Surface age','fontsize',8,'fontname','arial');
plot([0.4,0.7],[1040,1040],'k','linewidth',2)
text(0.9,1040,'S27 core','fontsize',8,'fontname','arial');
ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
text(0.4,780,'4','fontsize',8,'fontname','Arial','fontweight','bold')
text(13.8,780,'4''','fontsize',8,'fontname','Arial','fontweight','bold')
ylim([700 2200]);
%text(dist_bigblack_track4-0.04,-50,'115 ka','fontsize',8,'fontname','Arial','rotation',30);
%text(dist_bit16_track4-0.04,-50,'205 ka','fontsize',8,'fontname','Arial','rotation',30);
pos = get(gca,'pos');
set(gca,'pos',[0.1358 0.1307 0.6486 0.7443]);
set(gca,'fontsize',8,'fontname','arial')
set(gca,'xtick',[0:4:16])
set(gca,'ticklength',[0.02,0.05])
set(gca,'xticklabel',[])
text(0.3,2130,'b','fontsize',10,'fontname','arial','fontweight','bold');
export_fig(fullfile(REPO_HOME,'figures/track4_overview.pdf'),'-painters','-r600');
close;

% Plot age of picked layers
figure('units','inches');
hold off;
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 4.8 2.3],'color','white')
% Plot radar track4
uimagesc(distdata/1e3,zdata,filtdata,[-0.0015,0.0015]); hold on;
colormap(flipud(colormap('bone')));
set(gca,'ydir','normal');
freezeColors('nancolor',[1 1 1]);
hold on;
% Plot layer and bed picks
cm = cubehelix(20,1.0,-1.5,1,1,[0.2,0.8],[0.2,0.85]);
for i=1:size(track4.layages,1)
    if ~isnan(track4.layages(i))
      plot(track4.Distance_total4_all/1e3,track4.Height4_all'-track4.layheights(i,:),...
        'linewidth',1.5,'color',cm(round((track4.layages(i)-80)/10),:));
      %plot([dists_layers_outcropping_track4(i,1)/1e3,dists_layers_outcropping_track4(i,3)/1e3],...
      %   [dists_layers_outcropping_track4(i,2),0],'--','linewidth',2,'color',cm(round((track4.layages(i)-80)/10),:));
    end
end
plot(track4.Distance_total4_all(1:length(track4.thick))/1e3,track4.Height4_all-track4.thick,'k:','linewidth',1.5);
xlabel('Distance (km)','fontsize',8,'fontname','Arial')
ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
text(0.4,780,'4','fontsize',8,'fontname','Arial','fontweight','bold')
text(13.8,780,'4''','fontsize',8,'fontname','Arial','fontweight','bold')
ylim([700 2200]);
colormap(gca,cm); set(gca,'clim',[80 260]); h = colorbar; 
set(get(h,'title'),'string','Age (ka)','fontsize',8,'fontname','arial')
pos = get(get(h,'title'),'position');
set(get(h,'title'),'position',[pos(1)+1 pos(2) pos(3)])
set(h,'ticklength',[0.02 0.05])
pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2)-0.05 pos(3)+0.03 pos(4)]);
pos = get(h,'pos');
set(h,'pos',[pos(1) pos(2)-0.04 pos(3) pos(4)-0.01]);
set(gca,'fontsize',8,'fontname','arial')
set(gca,'xtick',[0:4:16])
set(gca,'ticklength',[0.02,0.05])
text(0.3,2130,'c','fontsize',10,'fontname','arial','fontweight','bold');
export_fig(fullfile(REPO_HOME,'figures/track4_ages.pdf'),'-painters','-r600');
close;

clear bedpicks h i pos thick ind laypicks


%% Interpolate ages from Track 4 to across tracks

 
%% Interpolate Track A
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_a_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_a_layers_max.splayer'));
thick = interp1(1:length(crosstracks.z_out),crosstracks.z_out,bedpicks)';
bed = [bed; crosstracks.Longitude1a(1:length(bedpicks)), crosstracks.Latitude1a(1:length(bedpicks)),...
       crosstracks.X1a(1:length(bedpicks)), crosstracks.Y1a(1:length(bedpicks)),...
       crosstracks.Height1a(1:length(bedpicks)),thick];
trackA.thick = thick;
labels_xy = [labels_xy;
             crosstracks.X1a(1)-800, crosstracks.Y1a(1)+100; 
             crosstracks.X1a(length(bedpicks))+150, crosstracks.Y1a(length(bedpicks))+200];
labels = [labels;
          'A ';'A'''];

% Get layer heights      
for i=1:size(laypicks,1)
      trackA.layheights(i,:) = interp1(trackA.z_out,laypicks(i,:))';
end

% Find where trackA intersects Track4
[~,ind_track4] = min(sqrt((track4.X4_all-trackA.interx_track4(1)).^2+(track4.Y4_all-trackA.interx_track4(2)).^2));
[~,ind_trackA] = min(sqrt((trackA.X-trackA.interx_track4(1)).^2+(trackA.Y-trackA.interx_track4(2)).^2));

% Get ages from Track 4
nonnan = find((~isnan(track4.layheights(:,ind_track4)) & (~isnan(track4.layages))));
trackA.layages_track4 = interp1(track4.layheights(nonnan,ind_track4),track4.layages(nonnan),trackA.layheights(:,ind_trackA));


%% Plot Track A 
plot_trackA = 0;
if plot_trackA
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(trackA.z_out==0,1,'last');
    uimagesc(trackA.Distance/1e3,trackA.z_out(ind:end),trackA.filtdata(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    % Plot layer and bed picks
    cm = cubehelix((260-80)/10,0.5,-1.5,1,1,[0.1,0.95],[0.2,0.9]);
    for i=1:size(trackA.layages_track4,1)
        if ~isnan(trackA.layages_track4(i))
            plot(trackA.Distance/1e3,trackA.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackA.layages_track4(i)-80)/10),:));
        %else
        %    plot(trackA.Distance/1e3,trackA.layheights(i,:),'color',[0.8 0 0],'linewidth',1);
        end
    end
    colormap(gca,cm); set(gca,'clim',[80 260]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    plot(trackA.Distance(1:length(trackA.thick))/1e3,trackA.thick,'k--','linewidth',2);
    % Plot intersecting transects
    %[~,ind] = min(sqrt((trackA.X-trackA.interx_track1B(1)).^1+(trackA.Y-trackA.interx_track1B(2)).^2));
    %plot([trackA.Distance(ind),trackA.Distance(ind)]/1e3,[0,trackA.z_out(end)],'k:','linewidth',1);
    %text(trackA.Distance(ind)/1e3-0.02,-50,'1','fontsize',8,'fontname','Arial');
    %[~,ind] = min(sqrt((trackA.X-trackA.interx_track2(1)).^2+(trackA.Y-trackA.interx_track2(2)).^2));
    %plot([trackA.Distance(ind),trackA.Distance(ind)]/1e3,[0,trackA.z_out(end)],'k:','linewidth',1);
    %text(trackA.Distance(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    %[~,ind] = min(sqrt((trackA.X-trackA.interx_track3(1)).^2+(trackA.Y-trackA.interx_track3(2)).^2));
    %plot([trackA.Distance(ind),trackA.Distance(ind)]/1e3,[0,trackA.z_out(end)],'k:','linewidth',1);
    %text(trackA.Distance(ind)/1e3-0.02,-50,'3','fontsize',8,'fontname','Arial');
    %[~,ind] = min(sqrt((trackA.X-trackA.interx_track4(1)).^2+(trackA.Y-trackA.interx_track4(2)).^2));
    %plot([trackA.Distance(ind),trackA.Distance(ind)]/1e3,[0,trackA.z_out(end)],'k:','linewidth',1);
    %text(trackA.Distance(ind)/1e3-0.02,-50,'4','fontsize',8,'fontname','Arial');
    text(0.05,1200,'A','fontsize',12,'fontname','Arial')
    text(2.5,1200,'A''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/trackA_picked.pdf'),'-painters','-transparent');
    close;
end
clear plot_trackA


%% Interpolate Track B
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_b_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_b_layers_max.splayer'));
thick = interp1(1:length(crosstracks.z_out),crosstracks.z_out,bedpicks)';
bed = [bed; crosstracks.Longitude1b(1:length(bedpicks)), crosstracks.Latitude1b(1:length(bedpicks)),...
       crosstracks.X1b(1:length(bedpicks)), crosstracks.Y1b(1:length(bedpicks)),...
       crosstracks.Height1b(1:length(bedpicks)),thick];
labels_xy = [labels_xy;
             crosstracks.X1b(1)+100, crosstracks.Y1b(1)+200; 
             crosstracks.X1b(length(bedpicks))-800, crosstracks.Y1b(length(bedpicks))];
labels = [labels;
          'B''';'B '];
trackB.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
      trackB.layheights(i,:) = interp1(trackB.z_out,laypicks(i,:))';
end

% Find where trackB intersects Track4
[~,ind_track4] = min(sqrt((track4.X4_all-trackB.interx_track4(1)).^2+(track4.Y4_all-trackB.interx_track4(2)).^2));
[~,ind_trackB] = min(sqrt((trackB.X-trackB.interx_track4(1)).^2+(trackB.Y-trackB.interx_track4(2)).^2));

% Get ages from Track 4
nonnan = find((~isnan(track4.layheights(:,ind_track4)) & (~isnan(track4.layages))));
trackB.layages_track4 = interp1(track4.layheights(nonnan,ind_track4),track4.layages(nonnan),trackB.layheights(:,ind_trackB));


%% Plot Track B
plot_trackB = 0;
if plot_trackB
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(trackB.z_out==0,1,'last');
    uimagesc(trackB.Distance/1e3,trackB.z_out(ind:end),trackB.filtdata(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(trackB.layheights,1)
        if ~isnan(trackB.layages_track4(i))
          plot(trackB.Distance/1e3,trackB.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackB.layages_track4(i)-80)/10),:));
        else
          plot(trackB.Distance/1e3,trackB.layheights(i,:),'r','linewidth',1);
        end
    end
    plot(trackB.Distance(1:length(trackB.thick))/1e3,trackB.thick,'k--','linewidth',2);
    [~,ind] = min(sqrt((trackB.X-trackB.interx_track1B(1)).^1+(trackB.Y-trackB.interx_track1B(2)).^2))
    plot([trackB.Distance(ind),trackB.Distance(ind)]/1e3,[0,trackB.z_out(end)],'k:','linewidth',1);
    text(trackB.Distance(ind)/1e3-0.02,-50,'1','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackB.X-trackB.interx_track2(1)).^2+(trackB.Y-trackB.interx_track2(2)).^2));
    plot([trackB.Distance(ind),trackB.Distance(ind)]/1e3,[0,trackB.z_out(end)],'k:','linewidth',1);
    text(trackB.Distance(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackB.X-trackB.interx_track3(1)).^2+(trackB.Y-trackB.interx_track3(2)).^2));
    plot([trackB.Distance(ind),trackB.Distance(ind)]/1e3,[0,trackB.z_out(end)],'k:','linewidth',1);
    text(trackB.Distance(ind)/1e3-0.02,-50,'3','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackB.X-trackB.interx_track4(1)).^2+(trackB.Y-trackB.interx_track4(2)).^2));
    plot([trackB.Distance(ind),trackB.Distance(ind)]/1e3,[0,trackB.z_out(end)],'k:','linewidth',1);
    text(trackB.Distance(ind)/1e3-0.02,-50,'4','fontsize',8,'fontname','Arial');
    text(0.05,1200,'B','fontsize',12,'fontname','Arial')
    text(2.5,1130,'B''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/trackB_picked.pdf'));
    close; 
end
clear plot_trackB


%% Interpolate Track C
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_c_bed.splayer'));
laypicks = load_splayer(fullfile('pickedlayers/layers_filtdata_c_layers_max.splayer'));
thick = interp1(1:length(trackC.z_out),trackC.z_out,bedpicks)';
bed = [bed; crosstracks.Longitude1c(1:length(bedpicks)), crosstracks.Latitude1c(1:length(bedpicks)),...
       crosstracks.X1c(1:length(bedpicks)), crosstracks.Y1c(1:length(bedpicks)),...
       crosstracks.Height1c(1:length(bedpicks)),thick];
labels_xy = [labels_xy;
             crosstracks.X1c(1)-800, crosstracks.Y1c(1); 
             crosstracks.X1c(length(bedpicks))+100, crosstracks.Y1c(length(bedpicks))+200];
labels = [labels;
          'C ';'C'''];
trackC.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
      trackC.layheights(i,:) = interp1(trackC.z_out,laypicks(i,:))';
end

% Find where trackC intersects Track4
[~,ind_track4] = min(sqrt((track4.X4_all-trackC.interx_track4(1)).^2+(track4.Y4_all-trackC.interx_track4(2)).^2));
[~,ind_trackC] = min(sqrt((trackC.X-trackC.interx_track4(1)).^2+(trackC.Y-trackC.interx_track4(2)).^2));

% Get ages from Track 4
nonnan = find((~isnan(track4.layheights(:,ind_track4)) & (~isnan(track4.layages))));
trackC.layages_track4 = interp1(track4.layheights(nonnan,ind_track4),track4.layages(nonnan),trackC.layheights(:,ind_trackC));


%% Plot Track C
plot_trackC = 0;
if plot_trackC
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(trackC.z_out==0,1,'last');
    uimagesc(trackC.Distance/1e3,trackC.z_out(ind:end),trackC.filtdata(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(trackC.layheights,1)
        if ~isnan(trackC.layages_track4(i))
          plot(trackC.Distance/1e3,trackC.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackC.layages_track4(i)-80)/10),:));
        else
          plot(trackC.Distance/1e3,trackC.layheights(i,:),'r','linewidth',1);
        end
    end
    plot(trackC.Distance(1:length(trackC.thick))/1e3,trackC.thick,'k--','linewidth',2);
    [~,ind] = min(sqrt((trackC.X-trackC.interx_track1B(1)).^1+(trackC.Y-trackC.interx_track1B(2)).^2));
    plot([trackC.Distance(ind),trackC.Distance(ind)]/1e3,[0,trackC.z_out(end)],'k:','linewidth',1);
    text(trackC.Distance(ind)/1e3-0.02,-50,'1','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackC.X-trackC.interx_track2(1)).^2+(trackC.Y-trackC.interx_track2(2)).^2));
    plot([trackC.Distance(ind),trackC.Distance(ind)]/1e3,[0,trackC.z_out(end)],'k:','linewidth',1);
    text(trackC.Distance(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackC.X-trackC.interx_track3(1)).^2+(trackC.Y-trackC.interx_track3(2)).^2));
    plot([trackC.Distance(ind),trackC.Distance(ind)]/1e3,[0,trackC.z_out(end)],'k:','linewidth',1);
    text(trackC.Distance(ind)/1e3-0.02,-50,'3','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackC.X-trackC.interx_track4(1)).^2+(trackC.Y-trackC.interx_track4(2)).^2));
    plot([trackC.Distance(ind),trackC.Distance(ind)]/1e3,[0,trackC.z_out(end)],'k:','linewidth',1);
    text(trackC.Distance(ind)/1e3-0.02,-50,'4','fontsize',8,'fontname','Arial');
    text(0.05,1200,'C','fontsize',12,'fontname','Arial')
    text(2.5,1130,'C''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/trackC_picked.pdf'));
    close;
end
clear plot_trackC

 
%% Interpolate Track D
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_d_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_d_layers_max.splayer'));
thick = interp1(1:length(crosstracks.z_out),crosstracks.z_out,bedpicks)';
bed = [bed; crosstracks.Longitude1d(1:length(bedpicks)), crosstracks.Latitude1d(1:length(bedpicks)),...
       crosstracks.X1d(1:length(bedpicks)), crosstracks.Y1d(1:length(bedpicks)),... 
       crosstracks.Height1d(1:length(bedpicks)),thick];
labels_xy = [labels_xy;
             crosstracks.X1d(1)+100, crosstracks.Y1d(1)+200; 
             crosstracks.X1d(length(bedpicks))-800, crosstracks.Y1d(length(bedpicks))];
labels = [labels;
          'D''';'D '];
trackD.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
    trackD.layheights(i,:) = interp1(trackD.z_out,laypicks(i,:))';
end

% Find where trackC intersects Track4
[~,ind_track4] = min(sqrt((track4.X4_all-trackD.interx_track4(1)).^2+(track4.Y4_all-trackD.interx_track4(2)).^2));
[~,ind_trackD] = min(sqrt((trackD.X-trackD.interx_track4(1)).^2+(trackD.Y-trackD.interx_track4(2)).^2));

% Get ages from Track 4
nonnan = find((~isnan(track4.layheights(:,ind_track4)) & (~isnan(track4.layages))));
trackD.layages_track4 = interp1(track4.layheights(nonnan,ind_track4),track4.layages(nonnan),trackD.layheights(:,ind_trackD));


%% Plot Track D
plot_trackD = 0;
if plot_trackD
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(crosstracks.z_out==0,1,'last');
    uimagesc(trackD.Distance/1e3,trackD.z_out(ind:end),trackD.filtdata(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(trackD.layheights,1)
        if ~isnan(trackD.layages_track4(i))
          plot(trackD.Distance/1e3,trackD.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackD.layages_track4(i)-80)/10),:));
        else
          plot(trackD.Distance/1e3,trackD.layheights(i,:),'r','linewidth',1);
        end
    end
    plot(trackD.Distance(1:length(trackD.thick))/1e3,trackD.thick,'k--','linewidth',2);
    [~,ind] = min(sqrt((trackD.X-trackD.interx_track1B(1)).^1+(trackD.Y-trackD.interx_track1B(2)).^2));
    plot([trackD.Distance(ind),trackD.Distance(ind)]/1e3,[0,trackD.z_out(end)],'k:','linewidth',1);
    text(trackD.Distance(ind)/1e3-0.02,-50,'1','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackD.X-trackD.interx_track2(1)).^2+(trackD.Y-trackD.interx_track2(2)).^2));
    plot([trackD.Distance(ind),trackD.Distance(ind)]/1e3,[0,trackD.z_out(end)],'k:','linewidth',1);
    text(trackD.Distance(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackD.X-trackD.interx_track3(1)).^2+(trackD.Y-trackD.interx_track3(2)).^2));
    plot([trackD.Distance(ind),trackD.Distance(ind)]/1e3,[0,trackD.z_out(end)],'k:','linewidth',1);
    text(trackD.Distance(ind)/1e3-0.02,-50,'3','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackD.X-trackD.interx_track4(1)).^2+(trackD.Y-trackD.interx_track4(2)).^2));
    plot([trackD.Distance(ind),trackD.Distance(ind)]/1e3,[0,trackD.z_out(end)],'k:','linewidth',1);
    text(trackD.Distance(ind)/1e3-0.02,-50,'4','fontsize',8,'fontname','Arial');
    text(2.4,1200,'D''','fontsize',12,'fontname','Arial');
    text(0.05,1200,'D','fontsize',12,'fontname','Arial');
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/trackD_picked.pdf'));
    close;    
end
clear plot_trackD
      
     
%% Interpolate Track E
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_e_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_e_layers_max.splayer'));
thick = interp1(1:length(crosstracks.z_out),crosstracks.z_out,bedpicks)';
bed = [bed; crosstracks.Longitude1e(1:length(bedpicks)), crosstracks.Latitude1e(1:length(bedpicks)),...
       crosstracks.X1e(1:length(bedpicks)), crosstracks.Y1e(1:length(bedpicks)),...
       crosstracks.Height1e(1:length(bedpicks)),thick];
labels_xy = [labels_xy;
             crosstracks.X1e(1)-800, crosstracks.Y1e(1); 
             crosstracks.X1e(length(bedpicks))+100, crosstracks.Y1e(length(bedpicks))+100];
labels = [labels;
          'E ';'E'''];
trackE.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
    trackE.layheights(i,:) = interp1(trackE.z_out,laypicks(i,:))';
end

% Find where trackC intersects Track4
trackE.interx_track4=interx([trackE.X,trackE.Y]',[track4.X4_all,track4.Y4_all]');
[~,ind_track4] = min(sqrt((track4.X4_all-trackE.interx_track4(1)).^2+(track4.Y4_all-trackE.interx_track4(2)).^2));
[~,ind_trackE] = min(sqrt((trackE.X-trackE.interx_track4(1)).^2+(trackE.Y-trackE.interx_track4(2)).^2));

% Get ages from Track 4
nonnan = find((~isnan(track4.layheights(:,ind_track4)) & (~isnan(track4.layages))));
trackE.layages_track4 = interp1(track4.layheights(nonnan,ind_track4),track4.layages(nonnan),trackE.layheights(:,ind_trackE));


%% Plot Track E
plot_trackE = 0;
if plot_trackE
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(trackE.z_out==0,1,'last');
    uimagesc(trackE.Distance/1e3,trackE.z_out(ind:end),trackE.filtdata(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(trackE.layheights,1)
        if ~isnan(trackE.layages_track4(i))
          plot(trackE.Distance/1e3,trackE.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackE.layages_track4(i)-80)/10),:));
        else
          plot(trackE.Distance/1e3,trackE.layheights(i,:),'r','linewidth',1);
        end
    end
    plot(trackE.Distance(1:length(trackE.thick))/1e3,trackE.thick,'k--','linewidth',2);
    [~,ind] = min(sqrt((trackE.X-trackE.interx_track1B(1)).^1+(trackE.Y-trackE.interx_track1B(2)).^2));
    plot([trackE.Distance(ind),trackE.Distance(ind)]/1e3,[0,trackE.z_out(end)],'k:','linewidth',1);
    text(trackE.Distance(ind)/1e3-0.02,-50,'1','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackE.X-trackE.interx_track2(1)).^2+(trackE.Y-trackE.interx_track2(2)).^2));
    plot([trackE.Distance(ind),trackE.Distance(ind)]/1e3,[0,trackE.z_out(end)],'k:','linewidth',1);
    text(trackE.Distance(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackE.X-trackE.interx_track3(1)).^2+(trackE.Y-trackE.interx_track3(2)).^2));
    plot([trackE.Distance(ind),trackE.Distance(ind)]/1e3,[0,trackE.z_out(end)],'k:','linewidth',1);
    text(trackE.Distance(ind)/1e3-0.02,-50,'3','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackE.X-trackE.interx_track4(1)).^2+(trackE.Y-trackE.interx_track4(2)).^2));
    plot([trackE.Distance(ind),trackE.Distance(ind)]/1e3,[0,trackE.z_out(end)],'k:','linewidth',1);
    text(trackE.Distance(ind)/1e3-0.02,-50,'4','fontsize',8,'fontname','Arial');
    text(0.05,1200,'E','fontsize',12,'fontname','Arial')
    text(2.3,1200,'E''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/trackE_picked.pdf'));
    close;       
end
clear plot_trackE

%%%%%%%%%%%%%%%


%% Interpolate Track 1B
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata1B_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_1b_layers_max.splayer'));
thick = interp1(1:length(track1B.z_out),track1B.z_out,bedpicks)';
bed = [bed; track1B.Longitude1B_all(1:length(bedpicks)), track1B.Latitude1B_all(1:length(bedpicks)),...
       track1B.X1B_all(1:length(bedpicks)), track1B.Y1B_all(1:length(bedpicks)),...
       track1B.Height1B_all(1:length(bedpicks)),thick];
labels_xy = [labels_xy; track1B.X1B_all(20)+200, track1B.Y1B_all(20)-500; 
             track1B.X1B_all(length(bedpicks)), track1B.Y1B_all(length(bedpicks))+500];
labels = [labels; '1 ';'1'''];
track1B.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
    track1B.layheights(i,:) = interp1(track1B.z_out,laypicks(i,:))';
end

% Find location of big black along transect
bigblack_track1B_loc=interx(bigblack',[track1B.X1B_all,track1B.Y1B_all]');
[~,ind_bigblack_track1B] = min(sqrt((track1B.X1B_all-bigblack_track1B_loc(1)).^2+(track1B.Y1B_all-bigblack_track1B_loc(2)).^2));
dist_bigblack_track1B = track1B.Distance_total1B_all(ind_bigblack_track1B)/1e3;

% Guess location where layers outcrop on the surface so that we can figure
% out which one is big black
dists_layers_outcropping_track1B = extrapolate_layers_to_surface(track1B.Distance_total1B_all,track1B.layheights);
% Find location of big black along transect
bigblack_track1B_loc=interx(bigblack',[track1B.X1B_all,track1B.Y1B_all]');
[~,ind_bigblack_track1B] = min(sqrt((track1B.X1B_all-bigblack_track1B_loc(1)).^2+(track1B.Y1B_all-bigblack_track1B_loc(2)).^2));
track1B.bigblack_dist = track1B.Distance_total1B_all(ind_bigblack_track1B);

% Get ages from transects A,B,C,D,E

% Find where track1B intersects trackA
[~,ind_track1B] = min(sqrt((track1B.X1B_all-track1B.interx_trackA(1)).^2+(track1B.Y1B_all-track1B.interx_trackA(2)).^2));
[~,ind_trackA] = min(sqrt((trackA.X-trackA.interx_track1B(1)).^2+(trackA.Y-trackA.interx_track1B(2)).^2));
% Get ages from TrackA
nonnan = find((~isnan(trackA.layheights(:,ind_trackA)) & (~isnan(trackA.layages_track4))));
track1B.layages_trackA = interp1(trackA.layheights(nonnan,ind_trackA),trackA.layages_track4(nonnan),track1B.layheights(:,ind_track1B));

% Find where track1B intersects trackB
[~,ind_track1B] = min(sqrt((track1B.X1B_all-track1B.interx_trackB(1)).^2+(track1B.Y1B_all-track1B.interx_trackB(2)).^2));
[~,ind_trackB] = min(sqrt((trackB.X-trackB.interx_track1B(1)).^2+(trackB.Y-trackB.interx_track1B(2)).^2));
% Get ages from TrackB
nonnan = find((~isnan(trackB.layheights(:,ind_trackB)) & (~isnan(trackB.layages_track4))));
track1B.layages_trackB = interp1(trackB.layheights(nonnan,ind_trackB),trackB.layages_track4(nonnan),track1B.layheights(:,ind_track1B));

% Find where track1B intersects trackC
[~,ind_track1B] = min(sqrt((track1B.X1B_all-track1B.interx_trackC(1)).^2+(track1B.Y1B_all-track1B.interx_trackC(2)).^2));
[~,ind_trackC] = min(sqrt((trackC.X-trackC.interx_track1B(1)).^2+(trackC.Y-trackC.interx_track1B(2)).^2));
% Get ages from TrackC
nonnan = find((~isnan(trackC.layheights(:,ind_trackC)) & (~isnan(trackC.layages_track4))));
track1B.layages_trackC = interp1(trackC.layheights(nonnan,ind_trackC),trackC.layages_track4(nonnan),track1B.layheights(:,ind_track1B));

% Find where track1B intersects trackD
[~,ind_track1B] = min(sqrt((track1B.X1B_all-track1B.interx_trackD(1)).^2+(track1B.Y1B_all-track1B.interx_trackD(2)).^2));
[~,ind_trackD] = min(sqrt((trackD.X-trackD.interx_track1B(1)).^2+(trackD.Y-trackD.interx_track1B(2)).^2));
% Get ages from TrackD
nonnan = find((~isnan(trackD.layheights(:,ind_trackD)) & (~isnan(trackD.layages_track4))));
track1B.layages_trackD = interp1(trackD.layheights(nonnan,ind_trackD),trackD.layages_track4(nonnan),track1B.layheights(:,ind_track1B));

% Find where track1B intersects trackE
[~,ind_track1B] = min(sqrt((track1B.X1B_all-track1B.interx_trackE(1)).^2+(track1B.Y1B_all-track1B.interx_trackE(2)).^2));
[~,ind_trackE] = min(sqrt((trackE.X-trackE.interx_track1B(1)).^2+(trackE.Y-trackE.interx_track1B(2)).^2));
% Get ages from TrackE
nonnan = find((~isnan(trackE.layheights(:,ind_trackE)) & (~isnan(trackE.layages_track4))));
track1B.layages_trackE = interp1(trackE.layheights(nonnan,ind_trackE),trackE.layages_track4(nonnan),track1B.layheights(:,ind_track1B));

% Set all ages
ages = [track1B.layages_trackA,track1B.layages_trackB,track1B.layages_trackC,track1B.layages_trackD,track1B.layages_trackE];
track1B.layages = nanmean(ages,2);
for i=1:size(track1B.layheights,1)
    nonnan = find(~isnan(ages(i,:)));
    if length(nonnan) > 1
        track1B.layerrors(i) = nanstd(ages(i,:))*1.96;
        track1B.laynums(i) = length(nonnan);
    else
        track1B.layerrors(i) = NaN;
        track1B.laynums(i) = NaN;
    end
end


%% Plot Track 1B
figure('units','inches');
hold off;
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 7 3],'color','white')
ind = find(track1B.z_out==0,1,'last');
uimagesc(track1B.Distance_total1B_all/1e3,track1B.z_out(ind:end),track1B.filtdata1B_all(ind:end,:),[-0.0015,0.0015]);
colormap(flipud(colormap('bone')));
xlabel('Distance (km)','fontsize',8,'fontname','Arial')
ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
freezeColors;
hold on;
for i=1:size(track1B.layages,1)
    if ~isnan(track1B.layages(i))
      plot(track1B.Distance_total1B_all/1e3,track1B.layheights(i,:),...
        'linewidth',1.5,'color',cm(round((track1B.layages(i)-80)/10),:));
    %else
      %[~,ind] = min(abs(track1B.Distance_total1B_all-14e3));
      %if track4.layheights(i,ind) < 600
      %  plot(track1B.Distance_total1B_all/1e3,track1B.layheights(i,:),'color',[0.8 0 0],'linewidth',1);
    end
end
colormap(gca,cm); set(gca,'clim',[80 260]); h = colorbar; 
set(get(h,'title'),'string','Age (ka)')
plot(track1B.Distance_total1B_all(1:length(track1B.thick))/1e3,track1B.thick,'k--','linewidth',2);
plot([dist_bigblack_track1B,dist_bigblack_track1B],[0,50],'w','linewidth',2);
text(dist_bigblack_track1B-0.04,-30,'115 ka','fontsize',8,'fontname','Arial','rotation',30);
%[~,ind] = min(sqrt((track1B.X1B_all-trackA.interx_track1B(1)).^2+(track1B.Y1B_all-trackA.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'A','fontsize',8,'fontname','Arial');
%[~,ind] = min(sqrt((track1B.X1B_all-trackB.interx_track1B(1)).^2+(track1B.Y1B_all-trackB.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'B','fontsize',8,'fontname','Arial');
%[~,ind] = min(sqrt((track1B.X1B_all-trackC.interx_track1B(1)).^2+(track1B.Y1B_all-trackC.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'C','fontsize',8,'fontname','Arial');
%[~,ind] = min(sqrt((track1B.X1B_all-trackD.interx_track1B(1)).^2+(track1B.Y1B_all-trackD.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'D','fontsize',8,'fontname','Arial');
%[~,ind] = min(sqrt((track1B.X1B_all-trackE.interx_track1B(1)).^2+(track1B.Y1B_all-trackE.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'E','fontsize',8,'fontname','Arial');
%[~,ind] = min(sqrt((track1B.X1B_all-trackF.interx_track1B(1)).^2+(track1B.Y1B_all-trackF.interx_track1B(2)).^2));
%plot([track1B.Distance_total1B_all(ind),track1B.Distance_total1B_all(ind)]/1e3,[0,track1B.z_out(end)],'k:','linewidth',1);
%text(track1B.Distance_total1B_all(ind)/1e3-0.02,-50,'F','fontsize',8,'fontname','Arial');
pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2)-0.03 pos(3) pos(4)-0.11]);
text(0.2,1150,'1','fontsize',12,'fontname','Arial')
text(19.2,1150,'1''','fontsize',12,'fontname','Arial')
colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
set(get(h,'title'),'string','Age (ka)')
ylim([0 1250]);
export_fig(fullfile(REPO_HOME,'figures/track1B_picked.pdf'));
close;  


%% Interpolate Track 2
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata2_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata2_layers_max.splayer'));
track2.thick = interp1(1:length(track2.z_out),track2.z_out,bedpicks)';
bed = [bed; track2.Longitude2_all(1:length(bedpicks)), track2.Latitude2_all(1:length(bedpicks)),...
       track2.X2_all(1:length(bedpicks)), track2.Y2_all(1:length(bedpicks)),...
       track2.Height2_all(1:length(bedpicks)),track2.thick];
labels_xy = [labels_xy;
             track2.X2_all(20)-200, track2.Y2_all(20)-500; 
             track2.X2_all(length(bedpicks))-650, -1.34468e6];
labels = [labels;
          '2 ';'2'''];

% Get layer heights      
for i=1:size(laypicks,1)
    track2.layheights(i,:) = interp1(track2.z_out,laypicks(i,:))';
end  

% Find location of big black along transect
bigblack_track2_loc=interx(bigblack',[track2.X2_all,track2.Y2_all]');
[~,ind_bigblack_track2] = min(sqrt((track2.X2_all-bigblack_track2_loc(1)).^2+(track2.Y2_all-bigblack_track2_loc(2)).^2));
dist_bigblack_track2 = track2.Distance_total2_all(ind_bigblack_track2)/1e3;

% Get ages from transects A,B,C,D,E

% Find where track2 intersects trackA
[~,ind_track2] = min(sqrt((track2.X2_all-track2.interx_trackA(1)).^2+(track2.Y2_all-track2.interx_trackA(2)).^2));
[~,ind_trackA] = min(sqrt((trackA.X-trackA.interx_track2(1)).^2+(trackA.Y-trackA.interx_track2(2)).^2));
% Get ages from TrackA
nonnan = find((~isnan(trackA.layheights(:,ind_trackA)) & (~isnan(trackA.layages_track4))));
track2.layages_trackA = interp1(trackA.layheights(nonnan,ind_trackA),trackA.layages_track4(nonnan),track2.layheights(:,ind_track2));

% Find where track2 intersects trackB
[~,ind_track2] = min(sqrt((track2.X2_all-track2.interx_trackB(1)).^2+(track2.Y2_all-track2.interx_trackB(2)).^2));
[~,ind_trackB] = min(sqrt((trackB.X-trackB.interx_track2(1)).^2+(trackB.Y-trackB.interx_track2(2)).^2));
% Get ages from TrackB
nonnan = find((~isnan(trackB.layheights(:,ind_trackB)) & (~isnan(trackB.layages_track4))));
track2.layages_trackB = interp1(trackB.layheights(nonnan,ind_trackB),trackB.layages_track4(nonnan),track2.layheights(:,ind_track2));

% Find where track2 intersects trackC
[~,ind_track2] = min(sqrt((track2.X2_all-track2.interx_trackC(1)).^2+(track2.Y2_all-track2.interx_trackC(2)).^2));
[~,ind_trackC] = min(sqrt((trackC.X-trackC.interx_track2(1)).^2+(trackC.Y-trackC.interx_track2(2)).^2));
% Get ages from TrackC
nonnan = find((~isnan(trackC.layheights(:,ind_trackC)) & (~isnan(trackC.layages_track4))));
track2.layages_trackC = interp1(trackC.layheights(nonnan,ind_trackC),trackC.layages_track4(nonnan),track2.layheights(:,ind_track2));

% Find where track2 intersects trackD
[~,ind_track2] = min(sqrt((track2.X2_all-track2.interx_trackD(1)).^2+(track2.Y2_all-track2.interx_trackD(2)).^2));
[~,ind_trackD] = min(sqrt((trackD.X-trackD.interx_track2(1)).^2+(trackD.Y-trackD.interx_track2(2)).^2));
% Get ages from TrackD
nonnan = find((~isnan(trackD.layheights(:,ind_trackD)) & (~isnan(trackD.layages_track4))));
track2.layages_trackD = interp1(trackD.layheights(nonnan,ind_trackD),trackD.layages_track4(nonnan),track2.layheights(:,ind_track2));

% Find where track2 intersects trackE
[~,ind_track2] = min(sqrt((track2.X2_all-track2.interx_trackE(1)).^2+(track2.Y2_all-track2.interx_trackE(2)).^2));
[~,ind_trackE] = min(sqrt((trackE.X-trackE.interx_track2(1)).^2+(trackE.Y-trackE.interx_track2(2)).^2));
% Get ages from TrackE
nonnan = find((~isnan(trackE.layheights(:,ind_trackE)) & (~isnan(trackE.layages_track4))));
track2.layages_trackE = interp1(trackE.layheights(nonnan,ind_trackE),trackE.layages_track4(nonnan),track2.layheights(:,ind_track2));

% Set all ages
ages = [track2.layages_trackA,track2.layages_trackB,track2.layages_trackC,track2.layages_trackD,track2.layages_trackE];
track2.layages = nanmean(ages,2);
for i=1:size(track2.layheights,1)
    nonnan = find(~isnan(ages(i,:)));
    if length(nonnan) > 1
        track2.layerrors(i) = nanstd(ages(i,:))*1.96;
        track2.laynums(i) = length(nonnan);
    else
        track2.layerrors(i) = NaN;
        track2.laynums(i) = NaN;
    end
end

      
%% Plot Track 2
plot_track2 = 0;
if plot_track2
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 7 3],'color','white')
    ind = find(track2.z_out==0,1,'last');
    uimagesc(track2.Distance_total2_all/1e3,track2.z_out(ind:end),track2.filtdata2_all(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(track2.layheights,1)
       if ~isnan(track2.layages(i))
          plot(track2.Distance_total2_all/1e3,track2.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((track2.layages(i)-80)/10),:));
       else
          plot(track2.Distance_total2_all/1e3,track2.layheights(i,:),'r','linewidth',1);
       end
    end
    plot(track2.Distance_total2_all(1:length(track2.thick))/1e3,track2.thick,'k--','linewidth',2);
    plot([dist_bigblack_track2,dist_bigblack_track2],[0,25],'color',[0 0 0.9],'linewidth',2);
    [~,ind] = min(sqrt((track2.X2_all-trackA.interx_track2(1)).^2+(track2.Y2_all-trackA.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'A','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track2.X2_all-trackB.interx_track2(1)).^2+(track2.Y2_all-trackB.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'B','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track2.X2_all-trackC.interx_track2(1)).^2+(track2.Y2_all-trackC.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'C','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track2.X2_all-trackD.interx_track2(1)).^2+(track2.Y2_all-trackD.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'D','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track2.X2_all-trackE.interx_track2(1)).^2+(track2.Y2_all-trackE.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'E','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track2.X2_all-trackF.interx_track2(1)).^2+(track2.Y2_all-trackF.interx_track2(2)).^2));
    plot([track2.Distance_total2_all(ind),track2.Distance_total2_all(ind)]/1e3,[0,track2.z_out(end)],'k:','linewidth',1);
    text(track2.Distance_total2_all(ind)/1e3-0.02,-50,'F','fontsize',8,'fontname','Arial');
    text(0.2,1200,'2','fontsize',12,'fontname','Arial')
    text(13.7,1200,'2''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/track2_picked.pdf'));
    close;   
end
clear plot_track2


%% Interpolate Track 3
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata3_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata3_layers_max.splayer'));
track3.thick = interp1(1:length(track3.z_out),track3.z_out,bedpicks)';
bed = [bed; track3.Longitude3_all(1:length(bedpicks)), track3.Latitude3_all(1:length(bedpicks)),...
       track3.X3_all(1:length(bedpicks)), track3.Y3_all(1:length(bedpicks)),...
       track3.Height3_all(1:length(bedpicks)),track3.thick];
labels_xy = [labels_xy;
             track3.X3_all(20)-400, track3.Y3_all(20)-500; 
             track3.X3_all(length(bedpicks))-1075, -1.34468e6];
labels = [labels;
          '3 ';'3'''];

% Get layer heights      
for i=1:size(laypicks,1)
    track3.layheights(i,:) = interp1(track3.z_out,laypicks(i,:))';
end         

% Find location of big black along transect
bigblack_track3_loc=interx(bigblack',[track3.X3_all,track3.Y3_all]');
[~,ind_bigblack_track3] = min(sqrt((track3.X3_all-bigblack_track3_loc(1)).^2+(track3.Y3_all-bigblack_track3_loc(2)).^2));
dist_bigblack_track3 = track3.Distance_total3_all(ind_bigblack_track3)/1e3;


% Get ages from transects A,B,C,D,E

% Find where track3 intersects trackA
[~,ind_track3] = min(sqrt((track3.X3_all-track3.interx_trackA(1)).^2+(track3.Y3_all-track3.interx_trackA(2)).^2));
[~,ind_trackA] = min(sqrt((trackA.X-trackA.interx_track3(1)).^2+(trackA.Y-trackA.interx_track3(2)).^2));
% Get ages from TrackA
nonnan = find((~isnan(trackA.layheights(:,ind_trackA)) & (~isnan(trackA.layages_track4))));
track3.layages_trackA = interp1(trackA.layheights(nonnan,ind_trackA),trackA.layages_track4(nonnan),track3.layheights(:,ind_track3));

% Find where track3 intersects trackB
[~,ind_track3] = min(sqrt((track3.X3_all-track3.interx_trackB(1)).^2+(track3.Y3_all-track3.interx_trackB(2)).^2));
[~,ind_trackB] = min(sqrt((trackB.X-trackB.interx_track3(1)).^2+(trackB.Y-trackB.interx_track3(2)).^2));
% Get ages from TrackB
nonnan = find((~isnan(trackB.layheights(:,ind_trackB)) & (~isnan(trackB.layages_track4))));
track3.layages_trackB = interp1(trackB.layheights(nonnan,ind_trackB),trackB.layages_track4(nonnan),track3.layheights(:,ind_track3));

% Find where track3 intersects trackC
[~,ind_track3] = min(sqrt((track3.X3_all-track3.interx_trackC(1)).^2+(track3.Y3_all-track3.interx_trackC(2)).^2));
[~,ind_trackC] = min(sqrt((trackC.X-trackC.interx_track3(1)).^2+(trackC.Y-trackC.interx_track3(2)).^2));
% Get ages from TrackC
nonnan = find((~isnan(trackC.layheights(:,ind_trackC)) & (~isnan(trackC.layages_track4))));
track3.layages_trackC = interp1(trackC.layheights(nonnan,ind_trackC),trackC.layages_track4(nonnan),track3.layheights(:,ind_track3));

% Find where track3 intersects trackD
[~,ind_track3] = min(sqrt((track3.X3_all-track3.interx_trackD(1)).^2+(track3.Y3_all-track3.interx_trackD(2)).^2));
[~,ind_trackD] = min(sqrt((trackD.X-trackD.interx_track3(1)).^2+(trackD.Y-trackD.interx_track3(2)).^2));
% Get ages from TrackD
nonnan = find((~isnan(trackD.layheights(:,ind_trackD)) & (~isnan(trackD.layages_track4))));
track3.layages_trackD = interp1(trackD.layheights(nonnan,ind_trackD),trackD.layages_track4(nonnan),track3.layheights(:,ind_track3));

% Find where track3 intersects trackE
[~,ind_track3] = min(sqrt((track3.X3_all-track3.interx_trackE(1)).^2+(track3.Y3_all-track3.interx_trackE(2)).^2));
[~,ind_trackE] = min(sqrt((trackE.X-trackE.interx_track3(1)).^2+(trackE.Y-trackE.interx_track3(2)).^2));
% Get ages from TrackE
nonnan = find((~isnan(trackE.layheights(:,ind_trackE)) & (~isnan(trackE.layages_track4))));
track3.layages_trackE = interp1(trackE.layheights(nonnan,ind_trackE),trackE.layages_track4(nonnan),track3.layheights(:,ind_track3));

% Set all ages
ages = [track3.layages_trackA,track3.layages_trackB,track3.layages_trackC,track3.layages_trackD,track3.layages_trackE];
track3.layages = nanmean(ages,2);
for i=1:size(track3.layheights,1)
    nonnan = find(~isnan(ages(i,:)));
    if length(nonnan) > 1
        track3.layerrors(i) = nanstd(ages(i,:))*1.96;
        track3.laynums(i) = length(nonnan);
    else
        track3.layerrors(i) = NaN;
        track3.laynums(i) = NaN;
    end
end


%% Plot Track 3
plot_track3 = 0;
if plot_track3
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 7 3],'color','white')
    ind = find(track3.z_out==0,1,'last');
    uimagesc(track3.Distance_total3_all/1e3,track3.z_out(ind:end),track3.filtdata3_all(ind:end,:),[-0.0015,0.0015]);
    colormap(flipud(colormap('bone')));
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    freezeColors;
    hold on;
    for i=1:size(track3.layheights,1)
       if ~isnan(track3.layages(i))
          plot(track3.Distance_total3_all/1e3,track3.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((track3.layages(i)-80)/10),:));
       else
          plot(track3.Distance_total3_all/1e3,track3.layheights(i,:),'r','linewidth',1);
       end
    end
    plot(track3.Distance_total3_all(1:length(track3.thick))/1e3,track3.thick,'k--','linewidth',2);
    plot([dist_bigblack_track3,dist_bigblack_track3],[0,25],'color',[0 0 0.9],'linewidth',2);
    [~,ind] = min(sqrt((track3.X3_all-trackA.interx_track3(1)).^2+(track3.Y3_all-trackA.interx_track3(2)).^2));
    plot([track3.Distance_total3_all(ind),track3.Distance_total3_all(ind)]/1e3,[0,track3.z_out(end)],'k:','linewidth',1);
    text(track3.Distance_total3_all(ind)/1e3-0.02,-50,'A','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track3.X3_all-trackB.interx_track3(1)).^2+(track3.Y3_all-trackB.interx_track3(2)).^2));
    plot([track3.Distance_total3_all(ind),track3.Distance_total3_all(ind)]/1e3,[0,track3.z_out(end)],'k:','linewidth',1);
    text(track3.Distance_total3_all(ind)/1e3-0.02,-50,'B','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track3.X3_all-trackC.interx_track3(1)).^2+(track3.Y3_all-trackC.interx_track3(2)).^2));
    plot([track3.Distance_total3_all(ind),track3.Distance_total3_all(ind)]/1e3,[0,track3.z_out(end)],'k:','linewidth',1);
    text(track3.Distance_total3_all(ind)/1e3-0.02,-50,'C','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track3.X3_all-trackD.interx_track3(1)).^2+(track3.Y3_all-trackD.interx_track3(2)).^2));
    plot([track3.Distance_total3_all(ind),track3.Distance_total3_all(ind)]/1e3,[0,track3.z_out(end)],'k:','linewidth',1);
    text(track3.Distance_total3_all(ind)/1e3-0.02,-50,'D','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((track3.X3_all-trackE.interx_track3(1)).^2+(track3.Y3_all-trackE.interx_track3(2)).^2));
    plot([track3.Distance_total3_all(ind),track3.Distance_total3_all(ind)]/1e3,[0,track3.z_out(end)],'k:','linewidth',1);
    text(track3.Distance_total3_all(ind)/1e3-0.02,-50,'E','fontsize',8,'fontname','Arial');
    text(0.2,1200,'3','fontsize',12,'fontname','Arial')
    text(14.7,1200,'3''','fontsize',12,'fontname','Arial')
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    set(get(h,'title'),'string','Age (ka)')
    export_fig(fullfile(REPO_HOME,'figures/track3_picked.pdf'));
    close;   
end
clear plot_track3
 

%% Interpolate Track F
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_1to2_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_1to2_layers_max.splayer'));
thick = interp1(1:length(trackF.z_out),trackF.z_out,bedpicks)';
bed = [bed; trackF.Longitude(1:length(bedpicks)), trackF.Latitude(1:length(bedpicks)),...
        trackF.X(1:length(bedpicks)), trackF.Y(1:length(bedpicks)),...
        trackF.Height(1:length(bedpicks)),thick];
trackF.thick = thick;
labels_xy = [labels_xy;
             trackF.X(20)-700, trackF.Y(20); 
             trackF.X(length(bedpicks))+50, trackF.Y(length(bedpicks))-400];
labels = [labels;
          'F ';'F'''];

% Get layer heights      
for i=1:size(laypicks,1)
    trackF.layheights(i,:) = interp1(trackF.z_out,laypicks(i,:))';
end

% Get ages for F from track1B
[~,ind_track1B] = min(sqrt((track1B.X1B_all-trackF.interx_track1B(1)).^2+(track1B.Y1B_all-trackF.interx_track1B(2)).^2));
[~,ind_trackF] = min(sqrt((trackF.X-trackF.interx_track1B(1)).^2+(trackF.Y-trackF.interx_track1B(2)).^2));
% Get ages from Track1B
nonnan = find((~isnan(track1B.layheights(:,ind_track1B)) & (~isnan(track1B.layages))));
trackF.layages_track1B = interp1(track1B.layheights(nonnan,ind_track1B),track1B.layages(nonnan),trackF.layheights(:,ind_trackF));

% Get ages for F from track2
[~,ind_track2] = min(sqrt((track2.X2_all-trackF.interx_track2(1)).^2+(track2.Y2_all-trackF.interx_track2(2)).^2));
[~,ind_trackF] = min(sqrt((trackF.X-trackF.interx_track2(1)).^2+(trackF.Y-trackF.interx_track2(2)).^2));
% Get ages from Track2
nonnan = find((~isnan(track2.layheights(:,ind_track2)) & (~isnan(track2.layages))));
trackF.layages_track2 = interp1(track2.layheights(nonnan,ind_track2),track2.layages(nonnan),trackF.layheights(:,ind_trackF));

ages = [trackF.layages_track1B,trackF.layages_track2];
trackF.layages = nanmean(ages,2);
for i=1:size(trackF.layheights,1)
    nonnan = find(~isnan(ages(i,:)));
    if length(nonnan) > 1
        trackF.layerrors(i) = nanstd(ages(i,:))*1.96;
        trackF.laynums(i) = length(nonnan);
    else
        trackF.layerrors(i) = NaN;
        trackF.laynums(i) = NaN;
    end
end


%% Plot Track F
plot_trackF = 0;
if plot_trackF
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 3],'color','white')
    ind = find(trackF.z_out==0,1,'last');
    uimagesc(trackF.Distance_total/1e3,trackF.z_out(ind:end),trackF.filtdata1a(ind:end,:),[-0.0015,0.0015]);
    hold on;
    colormap('bone');
    freezeColors;
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    for i=1:size(trackF.layheights,1)
        if ~isnan(trackF.layages(i))
            plot(trackF.Distance_total/1e3,trackF.layheights(i,:),...
            'linewidth',1.5,'color',cm(round((trackF.layages(i)-80)/10),:));
        else
            plot(trackF.Distance_total(1:length(trackF.layheights(i,:)))/1e3,trackF.layheights(i,:),'r','linewidth',1);
        end
    end
    plot(trackF.Distance_total(1:length(trackF.thick))/1e3,trackF.thick,'k--','linewidth',2);
    [~,ind] = min(sqrt((trackF.X-trackF.interx_track1B(1)).^2+(trackF.Y-trackF.interx_track1B(2)).^2));
    plot([trackF.Distance_total(ind),trackF.Distance_total(ind)]/1e3,[0,trackF.z_out(end)],'k:','linewidth',1);
    text(trackF.Distance_total(ind)/1e3-0.02,-50,'1B','fontsize',8,'fontname','Arial');
    [~,ind] = min(sqrt((trackF.X-trackF.interx_track2(1)).^2+(trackF.Y-trackF.interx_track2(2)).^2));
    plot([trackF.Distance_total(ind),trackF.Distance_total(ind)]/1e3,[0,trackF.z_out(end)],'k:','linewidth',1);
    text(trackF.Distance_total(ind)/1e3-0.02,-50,'2','fontsize',8,'fontname','Arial');
    ylim([0 1250]);
    colormap(gca,cm); set(gca,'clim',[80 255]); h = colorbar; 
    export_fig(fullfile(REPO_HOME,'figures/trackF_picked.pdf'));
    close;  
end
clear plot_trackF


%% Interpolate Track BIT58-1
BIT58_1 = load('BIT58_1.mat');
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_bit58_1_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_bit58_1_layers_max.splayer'));
thick = interp1(1:length(BIT58_1.z_out),BIT58_1.z_out,bedpicks)';
bed = [bed; BIT58_1.Longitude(1:length(bedpicks)), BIT58_1.Latitude(1:length(bedpicks)),...
        BIT58_1.X(1:length(bedpicks)), BIT58_1.Y(1:length(bedpicks)),...
        BIT58_1.Height(1:length(bedpicks)),thick];
BIT58_1.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
    BIT58_1.layheights(i,:) = interp1(BIT58_1.z_out,laypicks(i,:))';
end

[~,ind_BIT58_1] = min(sqrt((BIT58_1.X-labels_cores_xy(1,1)).^2+(BIT58_1.Y-labels_cores_xy(1,2)).^2));


%% Plot BIT58-1
plot_bit58_1 = 0;
if plot_bit58_1 
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 4 2.5],'color','white')
    ind = find(BIT58_1.z_out==0,1,'last');
    uimagesc(BIT58_1.Distance_total/1e3,BIT58_1.z_out(ind:end),BIT58_1.filtdata0(ind:end,:),[-0.0015,0.0015]);
    hold on;
    colormap('bone');
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    for i=1:size(BIT58_1.layheights,1)
       % if ~isnan(track1A.layages(i))
       %   plot(trackE.Distance/1e3,track1A.layheights(i,:),...
       %     'linewidth',1.5,'color',cm(round(track1A.layages(i))-80,:));
       % else
          plot(BIT58_1.Distance_total(1:length(BIT58_1.layheights(i,:)))/1e3,BIT58_1.layheights(i,:),'color',[0.8 0 0],'linewidth',1);
       % end
    end
    plot([BIT58_1.Distance_total(ind_BIT58_1),BIT58_1.Distance_total(ind_BIT58_1)]/1e3,[0,100],'w--','linewidth',3);
    plot(BIT58_1.Distance_total(1:length(BIT58_1.thick))/1e3,BIT58_1.thick,'k--','linewidth',2);
    ylim([0 1250]);
    export_fig(fullfile(REPO_HOME,'figures/bit58_1_picked.pdf'));
    close; 
end
clear plot_bit58_1


%% Interpolate Track BIT58-2
BIT58_2 = load('BIT58_2.mat');
bedpicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_bit58_2_bed.splayer'));
laypicks = load_splayer(fullfile(REPO_HOME,'pickedlayers/layers_filtdata_bit58_2_layers_max.splayer'));
thick = interp1(1:length(BIT58_2.z_out),BIT58_2.z_out,bedpicks)';
bed = [bed; BIT58_2.Longitude(1:length(bedpicks)), BIT58_2.Latitude(1:length(bedpicks)),...
        BIT58_2.X(1:length(bedpicks)), BIT58_2.Y(1:length(bedpicks)),...
        BIT58_2.Height(1:length(bedpicks)),thick];
BIT58_2.thick = thick;

% Get layer heights      
for i=1:size(laypicks,1)
    BIT58_2.layheights(i,:) = interp1(BIT58_2.z_out,laypicks(i,:))';
end

[~,ind_BIT58_2] = min(sqrt((BIT58_2.X-labels_cores_xy(1,1)).^2+(BIT58_2.Y-labels_cores_xy(1,2)).^2));


%% plot BIT58-2
plot_bit58_2 = 0;
if plot_bit58_2
    figure('units','inches');
    hold off;
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 6.5 2.5],'color','white')
    ind = find(BIT58_2.z_out==0,1,'last');
    uimagesc(BIT58_2.Distance_total/1e3,BIT58_2.z_out(ind:end),BIT58_2.filtdatadem(ind:end,:),[-0.0015,0.0015]);
    hold on;
    colormap('bone');
    xlabel('Distance (km)','fontsize',8,'fontname','Arial')
    ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    for i=1:size(BIT58_2.layheights,1)
       % if ~isnan(track1A.layages(i))
       %   plot(trackE.Distance/1e3,track1A.layheights(i,:),...
       %     'linewidth',1.5,'color',cm(round(track1A.layages(i))-80,:));
       % else
          plot(BIT58_2.Distance_total(1:length(BIT58_2.layheights(i,:)))/1e3,BIT58_2.layheights(i,:),'color',[0.8 0 0],'linewidth',1);
       % end
    end
    plot([BIT58_2.Distance_total(ind_BIT58_2),BIT58_2.Distance_total(ind_BIT58_2)]/1e3,[0,100],'w--','linewidth',3);
    %text(BIT58_2.Distance_total(ind_BIT58_2)+0.1,150,'BIT-58','fontsize',8,'fontname','arial','background','w');
    plot(BIT58_2.Distance_total(1:length(BIT58_2.thick))/1e3,BIT58_2.thick,'k--','linewidth',2);
    ylim([0 1250]);
    export_fig(fullfile(REPO_HOME,'figures/bit58_2_picked.pdf'),'-transparent',...
        '-painters','-p0.01');
    close;  
end
clear plot_bit58_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compile and save all bed picks 

% Data format: long lat x y zs thickness zb
bed = [bed,bed(:,5)-bed(:,6)];

% Save bed to text file
%dlmwrite('ahill_bed_picks.txt',[bed(:,1),bed(:,2),bed(:,5),bed(:,7)],'delimiter',' ','precision','%.6f');

clear thick bedpicks data


%% Make basemap

% Load landsat image
[I,xI,yI] = L8read('~/Desktop/AllanHills2016/Imagery/LC80591152016012LGN00','stretch');

figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 3.75 4.2],'color','white')
imagesc(xI,yI,I); hold on;
set(gca,'YDir','normal')
[x,y,Z] = bedmap2_data('surfacew',[5.104e5-11e3 5.175e5+3e3 5.175e5+3e3 5.104e5-11e3],[-1.3414e6+4e3 -1.3596e6-6e3 -1.3596e6-6e3 -1.3414e6+4e3],'xy');
[c,hc] = contour(x,y,Z,[1700:10:2200],'color',[0.5 0.5 0.5]);
%h(1)=plot(bigblack(:,1),bigblack(:,2),'-','color',[0.5 0.5 0.5],'linewidth',2);
%plot(bit16(:,1),bit16(:,2),'-','color',[0.5 0.5 0.5],'linewidth',2);
plot(stakes(:,1),stakes(:,2),'o','color','b','markersize',4);
quiver(stakes(:,1),stakes(:,2),stakes_u*5000,stakes_v*5000,0,'color','b','linewidth',1.5,'maxheadsize',5);
%plot(stakes(:,1),stakes(:,2),'ko','markersize',5,'markerfacecolor','b');
extend = dlmread('~/Dropbox/AllanHills2016/ShapeFiles/Flowlines/track1_extended.csv');
xtrack = [(extend(:,1))',flipud(track1B.X1B_all(8:end))']; % X-coord for flowline
ytrack = [(extend(:,2))',flipud(track1B.Y1B_all(8:end))']; % Y-coord for flowline
text(extend(1,1)-2.2e3,extend(1,2)+1.2e3,'Local Ice Divide','fontsize',8,'fontname','arial')
plot(xtrack,ytrack,'k:','linewidth',1.5)
plot(bed(:,3),bed(:,4),'k.','markersize',5)%,'markerfacecolor','flat')
%plot(xtrack,ytrack,'r-','linewidth',2);
for i=1:size(labels,1)
    %if strcmp(labels(i,:),'1''') || strcmp(labels(i,:),'1 ') || strcmp(labels(i,:),'4 ') ||...
    %        strcmp(labels(i,:),'4''')
        text(labels_xy(i,1),labels_xy(i,2),labels(i,:),'fontsize',8,'fontname','Arial');
    %end
end
for i=1:size(labels_cores,1)
    h(2) = plot(labels_cores_xy(i,1),labels_cores_xy(i,2),'ko','markerfacecolor','r','markersize',4);
end
text(labels_cores_xy(1,1)-2900,labels_cores_xy(1,2),labels_cores(1,:),'fontsize',8,'fontname','Arial');
text(labels_cores_xy(2,1)-2000,labels_cores_xy(2,2)+500,labels_cores(2,:),'fontsize',8,'fontname','Arial');
%plot(S27_surface_core_loc(:,1),S27_surface_core_loc(:,2),'kd','markerfacecolor','k');
axis('equal')
set(gca,'xtick',[5.0e5,5.05e5,5.10e5,5.15e5,5.2e5],'xticklabel',[500,505,510,515,520])
set(gca,'ytick',[-1365e3,-1360e3,-1355e3,-1350e3,-1345e3,-1340e3],'yticklabel',['','','','','',''])%[-1365,-1360,-1355,-1350,-1345,-1340])
text(499.9e3,-1361e3,'-1360','rotation',90,'fontsize',6,'fontname','arial')
text(499.9e3,-1356e3,'-1355','rotation',90,'fontsize',6,'fontname','arial')
text(499.9e3,-1351e3,'-1350','rotation',90,'fontsize',6,'fontname','arial')
text(499.9e3,-1346e3,'-1345','rotation',90,'fontsize',6,'fontname','arial')
text(499.9e3,-1341e3,'-1340','rotation',90,'fontsize',6,'fontname','arial')

xlabel('Easting (km)','fontsize',8,'fontname','Arial');
text(498.9e3,-1355e3,'Northing (km)','fontsize',8,'fontname','Arial','rotation',90)
set(gca,'fontsize',6,'fontname','Arial');
xlim([5.104e5-10e3,5.175e5+2e3]);
ylim([-1.3596e6-5e3,-1.3414e6+3e3])
%ylim([-1.3596e6-9e3,-1.3414e6+2e3])
pos = get(gca,'position');
set(gca,'position',[0.065 0.07 0.935 0.93])
box on
%Plot directional arrow towards south pole
[sp(1),sp(2)] = ll2ps(-90,0);
seed =[511000 -1343000];
sp_slope = (sp(2)-seed(2))/(sp(1)-seed(1));
quiver(seed(1),seed(2),-500,-500*sp_slope,0,'color','k','linewidth',2,'maxheadsize',5);
text(seed(1)-700,seed(2)-0.5e3,'South','fontsize',8,'fontname','Arial');
% Show zoom box for other figures
%plot([min(bed(:,3))-800,min(bed(:,3))-800,max(bed(:,3))+1000,max(bed(:,3))+1000,...
%    min(bed(:,3))-800],[min(bed(:,4))-1000,max(bed(:,4))+1000,max(bed(:,4))+1000,...
%    min(bed(:,4))-1000,min(bed(:,4))-1000],'k--');
% make scalebar
rectangle('Position',[500800,-1363200+14.5e3,8.0e3,5.5e3],'facecolor','w','edgecolor','k');
plot([501400,501400,503400,503400],[-1348400,-1348000,-1348000,-1348400],'k','linewidth',1.5);
plot([502400,502400],[-1348000,-1348250],'k','linewidth',1.5);
%plot(track4.X4_all(indices),track4.Y4_all(indices),'r','linewidth',1.5);
text(503900,-1348100,'2 km','fontsize',8,'fontname','Arial');
text(502000+1e3,-1359000+3e3,'Allan Hills','fontsize',8,'fontname','Arial','backgroundcolor','w');
quiver(501400,-1363200+19.5e3,0.2*5000,0,0,'color','b','linewidth',1.5,'maxheadsize',5)
text(502900,-1363200+19.5e3,'Velocity (0.2 m/a)','fontsize',8,'fontname','Arial');
%plot([501400,502500],[-1363200+19.5e3,-1363200+19.5e3],'-','color',[0.5 0.5 0.5],'linewidth',2)
%text(502900,-1363200+19.5e3,'Tephra layer','fontsize',8,'fontname','Arial');
plot([501400,502500],[-1363200+18.5e3,-1363200+18.5e3],'-','color','k','linewidth',2)
text(502900,-1363200+18.5e3,'Radar','fontsize',8,'fontname','Arial');
plot(501900,-1363200+17.5e3,'ko','markerfacecolor','r','markersize',4);
text(502900,-1363200+17.5e3,'Deep core','fontsize',8,'fontname','Arial');
plot(501900,-1363200+16.5e3,'o','color','b','markersize',4);
text(502900,-1363200+16.5e3,'Stake','fontsize',8,'fontname','Arial');
text(501000,-1363200+23.5e3,'a','fontsize',10,'fontname','Arial','fontweight','bold');

mapzoomps('Allan Hills','sw','insetsize',0.4);

% Add contour labels
%text(514.8e3,-1342.4e3,'2100 m','rotation',65,'fontsize',7,'fontname','arial')
%tl = clabel(c,hc,[1900:50:2200],'fontsize',7,'fontname','arial','labelspacing',500);
%tl = clabel(c,hc,'manual','fontsize',7,'fontname','arial','labelspacing',500);


export_fig(fullfile(REPO_HOME,'figures/ahills_basemap.pdf'),'-painters',...
    '-transparent','-r600');
close;

clear hc h sp seed sp_slope extend xtrack ytrack

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot bed elevations, ice thicknesses, and surface elevations

%dx = 50;
%xbed = 511000:dx:517600;
%ybed = -1359700:dx:-1345000;

%[xgrid,ygrid] = meshgrid(xbed,ybed);
%nonnan = find(~(isnan(bed(:,7))));
%[xcoord,ycoord] = ll2ps(bed(:,2),bed(:,1));
%zbed = griddata(xcoord(nonnan),ycoord(nonnan),bed(nonnan,7),xgrid,ygrid);
%for i=1:length(xbed)
%    for j=1:length(ybed)
%        [~,ind] = min(abs((xcoord-xbed(i)).^2+(ycoord-ybed(j)).^2));
%        if sqrt((xcoord(ind)-xbed(i)).^2+(ycoord(ind)-ybed(j)).^2) > 500
%            zbed(j,i) = NaN;
%        end
%    end
%end
%clear xcoord ycoord

% Bed elevation
fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 1.9 4.1],'color','white');
plot(0,0); hold on;

[x,y,Z] = bedmap2_data('bedw',[500000 520000 520000 500000],[-1300000 -1400000 -1400000 -1300000],'xy');
[c,ch] = contour(x,y,Z,[1000:100:2000],'k');
% label contours
text(516328,-1344735,'1000','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-90);
text(516362,-1342979,'1100','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-40);
text(513352,-1344986,'1200','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',75);
text(512599,-1346755,'1300','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',90);
text(512214,-1344255,'1400','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',80);
text(511262,-1346755,'1500','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',90);
text(510676,-1344255,'1600','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',80);
text(510157,-1353548,'1700','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',85);

scatter(bed(:,3),bed(:,4),5,bed(:,7)); hold on;
set(gca,'clim',[800,2000]);
colormap(flipud(cubehelix([],0.6,-0.9,1.9,1.1,[0,0.8],[0,0.9])))
h = colorbar('southoutside');
xlabel(h,'Bed elevation (m)','fontsize',8,'fontname','Arial');
set(h,'ticklength',[0.12 0.12]);
%contour(xbed,ybed,zbed,'k');
%freezeColors;
%h = imagesc(xI,yI,I); hold on;
%uistack(h,'bottom')
axis equal;
xlim([min(bed(:,3))-800,max(bed(:,3))+1000]);
ylim([min(bed(:,4))-1000,max(bed(:,4))+1000]);
set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
%set(fig,'Renderer','zbuffer');
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)-0.08 pos(3) pos(4)*0.5]);
set(h,'xtick',[1000 1400 1800]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1)+0.005 pos(2)+0.03 pos(3)+0.005 pos(4)]);
text(510200,-1341400,'c','fontweight','bold','fontname','Arial','fontsize',8);
% make scalebar
plot([514500,514500,516500,516500],[-1341200,-1340900,-1340900,-1341200],'k','linewidth',1.5)
text(514400,-1341700,'0','fontsize',8,'fontname','Arial');
text(516200,-1341700,'2 km','fontsize',8,'fontname','Arial');
for i=1:size(labels_cores,1)
    h(2) = plot(labels_cores_xy(i,1),labels_cores_xy(i,2),'ko','markerfacecolor','r','markersize',7);
end
text(labels_cores_xy(1,1)+500,labels_cores_xy(1,2)+200,labels_cores(1,:),'fontsize',8,'fontname','Arial');
text(labels_cores_xy(2,1)-1800,labels_cores_xy(2,2)+200,labels_cores(2,:),'fontsize',8,'fontname','Arial');

export_fig(fullfile(REPO_HOME,'figures/ahills_bed.pdf'),'-painters','-transparent')
close;

% Surface elevation
fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 1.9 4.1],'color','white');
plot(0,0); hold on;

[x,y,Z] = bedmap2_data('surfacew',[500000 520000 520000 500000],[-1300000 -1400000 -1400000 -1300000],'xy');
[c,hc] = contour(x,y,Z,[1800:25:2200],'k');
% label contours
text(511596,-1344066,'2100','fontsize',8,'fontname','Arial','backgroundcolor','w');
text(512116,-1348522,'2075','fontsize',8,'fontname','Arial','backgroundcolor','w');
text(510843,-1350505,'2050','fontsize',8,'fontname','Arial','backgroundcolor','w');
text(511262,-1351425,'2025','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-35);
text(510174,-1352044,'2000','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-30);
text(516780,-1357446,'1975','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',30);
text(513352,-1358533,'1950','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-40);
text(511354,-1359536,'1925','fontsize',8,'fontname','Arial','backgroundcolor','w');
text(510204,-1357505,'1950','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-90);
text(516800,-1359871,'1950','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',45);
text(511062,-1342276,'2125','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',55);

scatter(bed(:,3),bed(:,4),5,bed(:,7)+bed(:,6)); hold on;
set(gca,'clim',[1850,2150])
colormap(flipud(cubehelix([],0.6,-0.9,1.9,1.1,[0,0.8],[0,0.9])))
h = colorbar('southoutside');
xlabel(h,'Surface elevation (m)','fontsize',8,'fontname','Arial');
set(h,'ticklength',[0.12 0.12]);
%contour(xbed,ybed,zbed,'k');
%freezeColors;
%h = imagesc(xI,yI,I); hold on;
%uistack(h,'bottom')
axis equal;
xlim([min(bed(:,3))-800,max(bed(:,3))+1000]);
ylim([min(bed(:,4))-1000,max(bed(:,4))+1000]);
set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
%set(fig,'Renderer','zbuffer');
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)-0.08 pos(3) pos(4)*0.5]);
%set(h,'xtick',[1000 1400 1800]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1)+0.005 pos(2)+0.03 pos(3)+0.005 pos(4)]);
text(510200,-1341400,'a','fontweight','bold','fontname','Arial','fontsize',8);
% make scalebar
plot([514500,514500,516500,516500],[-1341200,-1340900,-1340900,-1341200],'k','linewidth',1.5)
text(514400,-1341700,'0','fontsize',8,'fontname','Arial');
text(516200,-1341700,'2 km','fontsize',8,'fontname','Arial');
for i=1:size(labels_cores,1)
    h(2) = plot(labels_cores_xy(i,1),labels_cores_xy(i,2),'ko','markerfacecolor','r','markersize',7);
end
text(labels_cores_xy(1,1)+500,labels_cores_xy(1,2)+200,labels_cores(1,:),'fontsize',8,'fontname','Arial');
text(labels_cores_xy(2,1)-1800,labels_cores_xy(2,2)+200,labels_cores(2,:),'fontsize',8,'fontname','Arial');

export_fig(fullfile(REPO_HOME,'figures/ahills_surface.pdf'),'-painters','-transparent');
close;

clear c ch x y Z h fig pos

% Ice thickness
fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 1.9 4.1],'color','white');
plot(0,0); hold on;

[x,y,Z1] = bedmap2_data('surfacew',[490000 530000 530000 490000],[-1300000 -1400000 -1400000 -1300000],'xy');
[x,y,Z2] = bedmap2_data('bedw',[490000 530000 530000 490000],[-1300000 -1400000 -1400000 -1300000],'xy');
[c,ch] = contour(x,y,Z1-Z2,[0:100:1200],'k');
% contour labels
text(516195,-1343147,'1000','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-50);
text(516212,-1344904,'1100','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',-70);
text(517059,-1349552,'          ','fontsize',6,'fontname','Arial','backgroundcolor','w','rotation',-45);
text(517359,-1349202,'1100','fontsize',8,'fontname','Arial','rotation',-50);
text(513519,-1344819,'900','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',80);
text(512683,-1346826,'800','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',90);
text(512014,-1344819,'700','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',85);
text(511445,-1346826,'600','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',90);
text(510760,-1344819,'500','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',85);
text(510191,-1346826,'400','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',90);
text(510174,-1353599,'300','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',85);
text(510258,-1359620,'200','fontsize',8,'fontname','Arial','backgroundcolor','w','rotation',95);

scatter(bed(:,3),bed(:,4),5,bed(:,6)); hold on;
colormap(flipud(cubehelix([],0.6,-0.9,1.9,1.1,[0,0.8],[0,0.9])))
h = colorbar('southoutside');
xlabel(h,'Ice thickness (m)','fontsize',8,'fontname','Arial');
set(gca,'clim',[0,1300]);
set(h,'ticklength',[0.12 0.12]);
%freezeColors;
%h = imagesc(xI,yI,I); hold on;
%uistack(h,'bottom')
axis equal;
xlim([min(bed(:,3))-800,max(bed(:,3))+1000]);
ylim([min(bed(:,4))-1000,max(bed(:,4))+1000]);
set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
%set(fig,'Renderer','zbuffer');
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)-0.08 pos(3) pos(4)*0.5]);
set(h,'xtick',[200 700 1200]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1)+0.005 pos(2)+0.03 pos(3)+0.005 pos(4)]);
text(510200,-1341400,'b','fontweight','bold','fontname','Arial','fontsize',8);
% make scalebar
plot([514500,514500,516500,516500],[-1341200,-1340900,-1340900,-1341200],'k','linewidth',1.5)
text(514400,-1341700,'0','fontsize',8,'fontname','Arial');
text(516200,-1341700,'2 km','fontsize',8,'fontname','Arial');
for i=1:size(labels_cores,1)
    h(2) = plot(labels_cores_xy(i,1),labels_cores_xy(i,2),'ko','markerfacecolor','r','markersize',7);
end
text(labels_cores_xy(1,1)+500,labels_cores_xy(1,2)+200,labels_cores(1,:),'fontsize',8,'fontname','Arial');
text(labels_cores_xy(2,1)-1800,labels_cores_xy(2,2)+200,labels_cores(2,:),'fontsize',8,'fontname','Arial');
%clabel(c,ch);

export_fig(fullfile(REPO_HOME,'figures/ahills_thickness.pdf'),'-painters','-transparent');
close;

clear c ch x y Z1 Z2 h fig pos
clear h pos xgrid ygrid dx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plots to show location of various radar data

plot_location_map = 0;
if plot_location_map
    [I,xI,yI] = L8read('~/Desktop/AllanHills2016/Imagery/LC80591152016012LGN00','stretch');

    fig = figure('units','inches');
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 1.9 3.9],'color','white');
    imagesc(xI,yI,I); hold on;
    set(gca,'YDir','normal')
    plot(bed(:,3),bed(:,4),'k.','markersize',5);
    %plot(trackA.X,trackA.Y,'r.','markersize',5);
    plot(track4.X4_all(ind_trackA:ind_trackD),track4.Y4_all(ind_trackA:ind_trackD),'r.','markersize',5);
    axis equal;
    xlim([min(bed(:,3))-800,max(bed(:,3))+1000]);
    ylim([min(bed(:,4))-1000,max(bed(:,4))+1000]);
    set(gca,'xtick',[],'ytick',[],'xticklabel',[],'yticklabel',[]);
    h(1)=plot(bigblack(:,1),bigblack(:,2),'-','color',[0.5 0.5 0.5],'linewidth',2);
    plot(bit16(:,1),bit16(:,2),'-','color',[0.5 0.5 0.5],'linewidth',2);
    for i=1:size(labels,1)
        text(labels_xy(i,1),labels_xy(i,2),labels(i,:),'fontsize',8,'fontname','Arial');
    end
    plot([514500,514500,516500,516500],[-1341200,-1340900,-1340900,-1341200],'k','linewidth',1.5)
    text(514400,-1341700,'0','fontsize',8,'fontname','Arial');
    text(516100,-1341700,'2 km','fontsize',8,'fontname','Arial');
    for i=1:size(labels_cores,1)
        plot(labels_cores_xy(i,1),labels_cores_xy(i,2),'ko','markerfacecolor','r','markersize',7);
    end
    text(labels_cores_xy(1,1)+500,labels_cores_xy(1,2)+200,labels_cores(1,:),'fontsize',8,'fontname','Arial');
    text(labels_cores_xy(2,1)-1800,labels_cores_xy(2,2)+200,labels_cores(2,:),'fontsize',8,'fontname','Arial');

    export_fig(fullfile(REPO_HOME,'figures/ahills_map.pdf'),'-painters','-transparent');
    close;
end

clear xI yI I pos fig plot_location_map


%% Determine age uncertainty

% First, we'll try just figuring out the average from all of the tracks
% where we have multiple age constraints. 

layages = [track1B.layages', track2.layages', track3.layages', trackF.layages'];
laynums = [track1B.laynums, track2.laynums, track3.laynums, trackF.laynums];
layerrors = [track1B.layerrors, track2.layerrors, track3.layerrors, trackF.layerrors];
nonnan = find(~isnan(layerrors));
laynums = laynums(nonnan);
layerrors = layerrors(nonnan);
layages = layages(nonnan);

% Make plot
figure('units','inches');
hold off;
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 3 3],'color','white')
for i=1:length(layerrors)
    plot(layages(i),layerrors(i),'ko','markerfacecolor','r','markersize',4); hold on;
end
xlabel('Age (ka)','fontsize',8,'fontname','arial');
ylabel('Average age difference (ka)','fontsize',8,'fontname','arial');
xlim([80,200])
ylim([0,7])
set(gca,'fontsize',8,'fontname','arial');
export_fig(fullfile(REPO_HOME,'figures/age_uncertainty.pdf'),'-painters',...
    '-transparent','-r600');
close;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preliminary Dansgaard Johnsen modeling (in ICE EQUIVALENT VALUES)
dansgaard_johnsen = 0;

if dansgaard_johnsen == 1
    
    % Possible drill locations
    [~,ind_trackA] = min(sqrt((track4.X4_all-trackA.interx_track4(1)).^2+(track4.Y4_all-trackA.interx_track4(2)).^2));
    [~,ind_trackD] = min(sqrt((track4.X4_all-trackD.interx_track4(1)).^2+(track4.Y4_all-trackD.interx_track4(2)).^2));
    [~,ind_1] = min(abs(track4.Distance_total4_all-14e3));
    ind_2 = ind_1+200;
    %indices = ind_core;%ind_trackA;%:ind_trackD;
    [~,indices] = min(abs(10-track4.Distance_total4_all/1e3));

    % For model sensitivity analysis
    bdots = 0.005:0.0001:0.015;
    hfracs = 0.1:0.01:0.7;

    misfits = zeros(length(bdots),length(hfracs),length(indices));
    depths_1ma = zeros(length(bdots),length(hfracs),length(indices));
    z = 0:1:max(track4.thick(indices));

    for j=1:length(indices)
        ind_core = indices(j);

        H = round(track4.thick(ind_core));
        ages = zeros(1,length(z));
        ages(:) = NaN;

        best_ages = zeros(length(z),length(hfracs));
        best_bdots = zeros(1,length(hfracs));
        best_depths_1ma = zeros(1,length(hfracs));
        best_depths_2ma = zeros(1,length(hfracs));
        best_ages(:,:) = NaN;
        best_bdots(:) = NaN;
        best_depths_1ma(:) = NaN;
        best_depths_2ma(:) = NaN;

        nonnan = find((~(isnan(track4.layheights(:,ind_core)))) & ...
                (~(isnan(track4.layages))));

        for l=1:length(hfracs);
            h = hfracs(l)*H;
            bdot_misfit_best=1e5;        

            for k=1:length(bdots)
                bdot = bdots(k);

                %Set values for H > z > h
                ind = find((z > h) & (z <= H));
                for i=ind
                    ages(i) = (2*H-h)/(2*bdot)*log((2*H-h)/(2*z(i)-h));  
                end
                [~,ind_h] = min(abs(z - h));
                ages(ind_h) = (2*H-h)/(2*bdot)*log((2*H-h)/h);
                ind = find(z < h);
                for i=ind
                    ages(i) = ages(ind_h)+(2*H-h)/bdot * (h/z(i)-1);
                end
                ages = ages/1000;

                % Find misfit
                modeled_ages = interp1(H-z,ages,track4.layheights(:,ind_core));
                misfits(k,l,j) = sum(((modeled_ages(nonnan)-track4.layages(nonnan))./track4.layages(nonnan)).^2);
                ind = find(~(isnan(ages)));
                depths_1ma(k,l,j) = interp1(ages(2:ind(end)),z(2:ind(end)),1000);

                if misfits(k,l,j) < bdot_misfit_best
                    bdot_misfit_best = misfits(k,l,j);
                    best_ages(:,l) = NaN;
                    best_ages(:,l) = ages;
                    best_bdots(l) = bdot;
                    best_depths_1ma(l) = interp1(ages(2:ind(end)),z(2:ind(end)),1000);
                    best_depths_2ma(l) = interp1(ages(2:ind(end)),z(2:ind(end)),2000);
                end

            end
        end

    end

    [~,ind1] = min(abs(hfracs - 0.5));
    [~,ind2] = min(abs(hfracs - 0.2));
    fig1 = figure('units','inches');
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 3 3],'color','white');
    h(2) = plot(best_ages(:,ind1),fliplr(z),'color',[0.7 0.7 0.7],'linewidth',2); hold on;
    h(1) = plot(best_ages(:,ind2),fliplr(z),'color',[0.2 0.2 0.2],'linewidth',2);
    plot(track4.layages,track4.layheights(:,ind_core),'.','color',[0 0 0.8],'markersize',10);
    set(gca,'ydir','reverse');
    ylim([0,H]);
    xlim([0 2.2e3])
    xlabel('Age (ka)','fontsize',8,'fontname','arial');
    ylabel('Ice thickness (m)','fontsize',8,'fontname','arial');
    pos1 = get(gca,'pos');
    set(gca,'ticklength',[0.020,0.05],'fontsize',8,'fontname','arial');
    [l1,l2] = legend(h,{'\lambda = 0.2','\lambda = 0.5'},'position',[0.33 0.8 0.1 0.08],'fontsize',8,'fontname','arial');
    plot([50,50,780,780,50],[20,190,190,20,20],'k');
    set(l2(3),'XData',[0.05 0.2])
    pos = get(l2(1),'position');
    set(l2(1),'position',[0.25 pos(2) pos(3)])
    set(l2(5),'XData',[0.05 0.2])
    pos = get(l2(2),'position');
    set(l2(2),'position',[0.25 pos(2) pos(3)])
    legend('boxoff');
    plot([60,60,225,225,60],[450,910,910,450,450],'k--');
    text(-300,H,'Bed','fontsize',8,'fontname','arial');

    fig2 = figure('units','inches');
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 0.225 0.877],'color','white');
    plot(best_ages(:,ind1),fliplr(z),'color',[0.7 0.7 0.7],'linewidth',2); hold on;
    plot(best_ages(:,ind2),fliplr(z),'color',[0.2 0.2 0.2],'linewidth',2);
    plot(track4.layages,track4.layheights(:,ind_core),'.','color',[0 0 0.8],'markersize',10);
    xlim([60 225])
    ylim([450 910])
    set(gca,'ydir','reverse');
    set(gca,'xtick',[100:100:500],'fontsize',8)
    set(gca,'ytick',[0:200:1000],'fontsize',8);
    set(gca,'ticklength',[0.040,0.08])

    inset(fig1,fig2,0.3,0.8);
    export_fig(fullfile(REPO_HOME,'figures/ahills_djmodel.pdf'),'-painters',...
        '-transparent');
    close;
    clear fig1 fig2 l2 ind1 ind2 h

    fig = figure('units','inches');
    pos = get(gcf,'pos');
    set(gcf,'pos',[0 0 3 3],'color','white');
    [c,h] = contour(hfracs,bdots*1e3,log(misfits(:,:,1)),'linewidth',1);
    colormap('bone')
    clabel(c,h,[-1,0,1,2,3]); hold on;
    filtlen = 3;
    best_bdots_filt = zeros(1,length(best_bdots));
    best_bdots_filt(1:filtlen) = best_bdots(1:filtlen);
    best_bdots_filt(end-filtlen:end) = best_bdots(end-filtlen:end);
    for i = filtlen+1:length(best_bdots)-filtlen-1
        best_bdots_filt(i) = mean(best_bdots(i-filtlen:i+filtlen));
    end

    plot(hfracs,best_bdots_filt*1e3,'color',[0.8 0 0],'linewidth',2);
    ylabel('Accumulation rate (mm/yr)','fontsize',8,'fontname','Arial');
    xlabel('\lambda','fontsize',8,'fontname','Arial');
    %set(gca,'fontsize',8,'fontname','Arial');
    export_fig(fullfile(REPO_HOME,'figures/ahills_djmodel_sensitivity.pdf'),...
        '-painters','-transparent');
    close;
    clear c h filtlen fig pos

end

clear nonnan laypicks

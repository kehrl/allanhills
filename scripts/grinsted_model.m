% This script runs the Grinsted et al 2003 flowline model for different
% combinations of accumulation rate "bdot_constants" and velocity decrease during
% glacial periods "us_constants." For each combination of variables, it 
% produces a set of figures. The model can use Radar Track 1 or 4 as the 
% flowline. 
%
% Laura Kehrl, UW, 10/01/2017

% Model inputs
%bdot_constant = 0.009; %mm/yr we
%us_constant = 0.25;

bdot_constants = [0.007];%:0.0005:0.011];
us_constants = 0.8;%[0.0:0.10:1.00];

enhanced_bdot = 0;
save_rmse_file = 0;

if save_rmse_file
    fid = fopen(fullfile(REPO_HOME,sprintf('RMSE/RMSE_bdot%02dmm.txt',bdot_constants(1)*1000),'w'));
end

for m=1:length(bdot_constants)
for n=1:length(us_constants)

%% Set up model
bdot_constant = bdot_constants(m);
us_constant = us_constants(n);

disp(sprintf('Running model with \n Us = %d%% \n bdot = %d mm',us_constant*100,bdot_constant*1000));

% Which flowline?
flowline = 'track1';

clear H_ice H_ice_temp

if strcmp(flowline,'track4')
    % Extended flowline
    extend = dlmread('~/Dropbox/AllanHills2016/ShapeFiles/Flowlines/flowline_extended.csv');
    extend = flipud(extend);

    % Inputs
    xtrack = [(extend(:,1))',flipud(track4.X4_all(13:end))']; % X-coord for flowline
    ytrack = [(extend(:,2))',flipud(track4.Y4_all(13:end))']; % Y-coord for flowline
    dist = distance_along_transect(xtrack,ytrack); % New dist coord for flowline
    dist_offset = max(dist) - max(track4.Distance_total4_all);
    thick_ice = [bedmap2_interp(extend(:,1),extend(:,2),'surfacew')'-bedmap2_interp(extend(:,1),...
        extend(:,2),'bedw')',flipud(track4.thick(13:end))']; % Ice thickness
    zs = [bedmap2_interp(extend(:,1),extend(:,2),'surfacew')',...
        flipud(track4.Height4_all(13:end))']; % Surface elevation (wgs 84)
elseif strcmp(flowline,'track1')
    extend = dlmread('~/Dropbox/AllanHills2016/ShapeFiles/Flowlines/track1_extended.csv');
    extend = extend;
    
    % Inputs
    xtrack = [(extend(:,1))',flipud(track1B.X1B_all(8:end))']; % X-coord for flowline
    ytrack = [(extend(:,2))',flipud(track1B.Y1B_all(8:end))']; % Y-coord for flowline
    dist = distance_along_transect(xtrack,ytrack); % New dist coord for flowline
    dist_offset = max(dist) - max(track1B.Distance_total1B_all);
    thick_ice = [bedmap2_interp(extend(:,1),extend(:,2),'surfacew')'-bedmap2_interp(extend(:,1),...
        extend(:,2),'bedw')',flipud(track1B.thick(8:end))']; % Ice thickness
    ind = find(~(isnan(thick_ice)));
    thick_ice = interp1(dist(ind),thick_ice(ind),dist)';
    zs = [bedmap2_interp(extend(:,1),extend(:,2),'surfacew')',...
        flipud(track1B.Height1B_all(8:end))']; % Surface elevation (wgs 84) 
        
end

% Constants
rho_i = 900; % Ice density
rho_w = 1000; % Water density

% Velocity shapefactor
%k = 5; % linear
%k = 3.4; % constant 
k = 5;

% Grid spacing
dx = 250;
dy = 0.001;

% Normalized vertical direction
zbar = 0:dy:1.0;

% Velocity profile shapefactor
f = tanh(k*zbar)/tanh(k);

% Model timesteps
dt = 500; % years
t = 0:dt:1.5e6; % time

% Model variables
x = [dist(1):dx:dist(end)]; % distance along grid
zs_x = interp1(dist,zs,x); % surface elevation along grid
H_ice = interp1(dist,thick_ice,x); % ice thickness along grid
H_we = H_ice*(rho_i/rho_w); % W.E. ice thickness along grid 

% We want to decrease velocities and accumulation rates during glacial
% periods. To do this, we use dates from Masson-Delmotte, et al. 
% EPICA DOME C record of glacial and interglacial intensities (2010). Note
% that it is assumed that accumulation rates decrease by 50% during glacial
% periods, which is similar to the decrease observed in the EPICA core (see
% below).
Us_t = ones(length(t),1);
bdot_t = ones(length(t),1);
for i = 1:length(t)
    if ((t(i)/1e3 > 11.6) && (t(i)/1e3 < 115.0)) || ...
            ((t(i)/1e3 > 132.0) && (t(i)/1e3 < 233.6)) || ...
            ((t(i)/1e3 > 245.6) && (t(i)/1e3 < 316.8)) || ...
            ((t(i)/1e3 > 335.6) && (t(i)/1e3 < 392.0)) || ...
            ((t(i)/1e3 > 424.6) && (t(i)/1e3 < 470.0)) || ...
            ((t(i)/1e3 > 529.6) && (t(i)/1e3 < 560.0)) || ...
            ((t(i)/1e3 > 580.0) && (t(i)/1e3 < 603.0)) || ...
            ((t(i)/1e3 > 628.0) && (t(i)/1e3 < 680.0)) || ...
            ((t(i)/1e3 > 720.0) && (t(i)/1e3 < 769.6)) || ...
            ((t(i)/1e3 > 791.6))
        bdot_t(i) = 0.5;
        Us_t(i) = us_constant;
    end
end
% Compare to EPICA Dome C accumulation rates
epica = xlsread('EPICA/AICC2012_official.xls','EDC');

ind_int = find(bdot_t==1);
ind_glac = find(bdot_t ~= 1);
epica_int_bdot = interp1(epica(:,2),epica(:,6),t(ind_int));
epica_glac_bdot = interp1(epica(:,2),epica(:,6),t(ind_glac));
% Decrease during glacials given by 
% nanmean(epica_glac_bdot)/nanmean(epica_int_bdot)
% i.e., accumulation rates were ~30% less during glacial periods


%% Get surface velocities and accumulation/ablation rates 
% Get values from stakes (A01, P04, P06, P15, P18)
n=1;
clear bdot_pts x_pts Us_pts
for i=1:length(stakes)
    [mindist,j] = min(sqrt((stakes(i,1)-xtrack).^2+(stakes(i,2)-ytrack).^2));
    % Get data from stakes within 0.5 km of flowline
    if mindist < 0.5e3
        bdot_pts(n) = stakes(i,10);
        Us_pts(n) = stakes(i,4);
        x_pts(n) = dist(j);
        n = n+1;
    end
end
bdot_pts = bdot_pts*rho_i/rho_w; % convert bdot to w.e.

% We don't have many velocity measurements, so one option is to interpolate
% between the observations using the shallow ice approximation. Note that
% I'm not presently using this approach, because it does not greatly reduce
% the RMS (so the simpler option of linear interpolation is better).
x_temp = [dist(1):2.5:dist(end)];
zs_temp = interp1(dist,zs,x_temp);
H_ice_temp = interp1(dist,thick_ice,x_temp);

% Smooth with a low pass filter
[b,a]=butter(4,0.001,'low');
zs_temp_mean = filtfilt(b,a,zs_temp);

dhdx = gradient(zs_temp_mean,x_temp);
d2hdx2 = gradient(dhdx,x_temp);
g = 9.8;
tau_d = -1*(rho_i*g*dhdx.*H_ice_temp);
ind = find((tau_d < 0));
ind2 = find((tau_d > 0));
tau_d(ind) = mean(tau_d(ind2));
clear ind ind2

% Get front multiplier from existing data, assuming a linear viscous
% material
if strcmp(flowline,'track1') || strcmp(flowline,'track4')
    taud_pts = abs(interp1(x_temp,tau_d,x_pts));
    H_pts = interp1(x_temp,H_ice_temp,x_pts);
    fracs = Us_pts./((taud_pts).*H_pts);
    frac = nanmean(fracs);
    Us_temp = frac*(tau_d).*H_ice_temp;
    Us_temp_mean = zeros(size(Us_temp));
    Us_temp_mean(:) = NaN;
    for i = 101:length(x_temp)-101
        Us_temp_mean(i) = mean(Us_temp(i-50:i+50));
    end
    Us_temp_mean(1) = 0.001;
end

% Adjut parameters
x_pts(end+1) = 0; % ice divide
bdot_pts(end+1) = bdot_constant; % for a constant accumulation rate upstream
Us_pts(end+1) = 0.000; % velocities decrease to near 0 at the ice divide

% Create timeseries of surface velocity Us
% Use Us proportional to tau_d
% ind = find(~isnan(Us_temp_mean));
% Us_back = interp1(x_temp(ind),Us_temp_mean(ind),x,'linear','extrap');
% Use linear velocity profile, using measured velocities
ind = find(~isnan(Us_pts));
Us_back = interp1([x_pts(ind),x(end)],[Us_pts(ind),Us_pts(ind(2))],x,'linear','extrap');
Us = zeros(length(x),length(t));
for j=1:length(t)
    Us(:,j) = Us_back*Us_t(j);
end

clear frac fracs dhdx H_pts taud_pts tau_d zs_temp_mean Us_temp_mean x_temp

% Get surface slope and curvature
x_temp = fliplr([track1B.Distance_total1B_all(1):2.5:track1B.Distance_total1B_all(end)]);
zs_temp = interp1(fliplr(track1B.Distance_total1B_all),track1B.Height1B_all,x_temp);

[b,a]=butter(4,0.0011,'low');
zs_temp(1:6) = zs_temp(7);
zs_temp_mean = filtfilt(b,a,zs_temp);

dhdx = gradient(zs_temp,x_temp);
dhdx_mean = filtfilt(b,a,dhdx);
d2hdx2_mean = (gradient(dhdx_mean,x_temp));

d2hdx2 = interp1(x_temp+dist_offset,d2hdx2_mean,x);

% Create timeseries of bdot
ind = find(~(isnan(bdot_pts)));
bdot_back = interp1([x_pts(ind),16e3,17e3,24e3],[bdot_pts(ind),bdot_constant,0,-0.02],x,'linear','extrap');
bdot = zeros(length(x),length(t));
if enhanced_bdot
    ind = find((x >= 6.0e3) & (x <= 6.5e3));
    bdot_back(ind) = bdot_back(ind)*2;
end
ind = find((bdot_back > 0));
%bdot_back(ind) = bdot_back(ind).*(1+d2hdx2(ind)*1e5);
for j=1:length(t)
    bdot(:,j) = bdot_back;
    bdot(ind,j) = bdot_back(ind)*bdot_t(j);
end


%% Set up horizontal and vertical velocities

% Assuming model is along a flowline, so there should be no outward flow
Vs = 0;

% Model velocities;
U = zeros(length(zbar),length(x),length(t));
for i=1:length(x)
    for j = 1:length(t)
        U(:,i,j) = f*Us(i,j);
    end
end
V = 0;

% Horizontal velocity gradients
dUsdx = zeros(length(x),length(t));
for i = 2:length(x)-1
    for j = 1:length(t)
        dUsdx(i,j) = (Us(i+1,j)-Us(i-1,j))/(2*dx);
    end
end
dUsdx(1,:) = dUsdx(2,:);
dUsdx(end) = dUsdx(end-1);
dUdx = f'*dUsdx(:,1)';

dVsdy = zeros(length(x),length(t));
for j =1:length(t)
    dVsdy(:,j) = (bdot(:,j)./(mean(f)*H_we'))-dUsdx(:,j);
end
dVdy = f'*dVsdy(:,1)';

% Vertical velocity gradient
dWdz = zeros(length(zbar),length(x),length(t));
for i=1:length(x)
    for j=1:length(t)
        dWdz(:,i,j) = -(dUsdx(i,j)+dVsdy(i,j))*f;
    end
end
    
% Normalized vertical velocity    
wbar = zeros(length(zbar),length(x),length(t));
for i = 1:length(x)
    for j = 1:length(zbar)
        for n =1:length(t)
            wbar(j,i,n) = ((dUsdx(i,n) + dVsdy(i,n))/(-k*tanh(k))) * (log((1+exp(-2*k*zbar(j)))/2) + k*zbar(j));
        end
    end
end


%% Get ages at x=0 from Dansgaard Johnsen model
% Get ages at ice divide assuming a simple 1D Dansgaard Johnsen model. This
% serves as the upstream boundary condition if a point is tracked out of
% the upstream vertical boundary (the "ice divide").
H_we_divide = round(H_we(1));
h_we_divide = 0.7*H_we_divide;
z_we_divide = 0:0.1:H_we_divide;
ages_we_divide = zeros(1,length(z_we_divide));
ages_we_divide = nan;
bdot_we_divide = bdot(1);

ind = find(z_we_divide > h_we_divide);
for i=ind
	ages_we_divide(i) = (2*H_we_divide-h_we_divide)/(2*bdot_we_divide)*log((2*H_we_divide-h_we_divide)/(2*z_we_divide(i)-h_we_divide));  
end
[~,ind_h] = min(abs(z_we_divide - h_we_divide));
ages_we_divide(ind_h) = (2*H_we_divide-h_we_divide)/(2*bdot_we_divide)*log((2*H_we_divide-h_we_divide)/h_we_divide);
[~,ind] = find(z_we_divide < h_we_divide);
for i=ind
    ages_we_divide(i) = ages_we_divide(ind_h)+(2*H_we_divide-h_we_divide)/bdot_we_divide * (h_we_divide/z_we_divide(i)-1);
end


%% Lagrangian integration to get age for each grid cell

ages = zeros(length(zbar),length(x));
%ages(:,1) = interp1(z_we_divide,ages_we_divide,zbar*H_we(1));
for j = 1:length(zbar)
    for n = 1:length(x)
        pt(1,1) = x(n);
        pt(1,2) = zbar(j);
        pt(1,3) = interp1(x,zs_x,pt(1,1))-interp1(x,H_we,pt(1,1))*rho_w/rho_i+pt(1,2)*rho_w/rho_i*interp1(x,H_we,pt(1,1));
        pointwithindomain = 1;

        while pointwithindomain
            for i = 1:length(t)
                u_pt_0 = interp2(x,zbar,U(:,:,i),pt(i,1),pt(i,2),'cubic');
                v_pt_0 = interp2(x,zbar,wbar(:,:,i),pt(i,1),pt(i,2),'cubic');
                xpt_1 = pt(i,1)-u_pt_0*dt;
                zpt_1 = pt(i,2)-v_pt_0*dt;
                H_we_pt_1 = interp1(x,H_we,xpt_1);
                zs_pt_1 = interp1(x,zs_x,xpt_1);
                if (xpt_1 <= max(x)) & (xpt_1 >= 0) & (zpt_1 <= 1) & (zpt_1 >= 0) 
                    pt(i+1,1) = xpt_1;
                    pt(i+1,2) = zpt_1;
                    pt(i+1,3) = zs_pt_1-H_we_pt_1*rho_w/rho_i + zpt_1*rho_w/rho_i*H_we_pt_1;
                else
                    if xpt_1 <= 0 & i == 2
                        ages(j,n) = t(i);
                    elseif xpt_1 <= 0
                        ages(j,n) = t(i);
                    end
                    pointwithindomain = 0;
                    % Quit while loop if the point leaves the glacier
                    % domain
                    break
                
                end
                clear zs_pt_1 H_we_pt_1 xpt_1 zpt_1 u_pt_0 v_pt_0  
            end
            pointwithindomain = 0;
        end
        clear pointwithindomain
        if i ~=2
            % If the point moves back in time, calculate how old it is 
            ages(j,n) = t(i); 
        end
        if (x(n)/1e3 == 17) && (j == length(zbar))
           pt17 = pt;
        end
        if (x(n)/1e3 == 18) && (j == length(zbar))
           pt18 = pt;
        end
        if (x(n)/1e3 == 19) && (j == length(zbar))
           pt19 = pt;
        end
        if (x(n)/1e3 == 20) && (j == length(zbar))
           pt20 = pt;
        end
        if (x(n)/1e3 == 21) && (j == length(zbar))
           pt21 = pt;
        end
        if (x(n)/1e3 == 22) && (j == length(zbar))
           pt22 = pt;
        end 
        if (x(n)/1e3 == 23) && (j == length(zbar))
           pt23 = pt;
        end
        if (x(n)/1e3 == 24) && (j == length(zbar))
           pt24 = pt;
        end
        clear pt
    end
end


%% Make new grid in x, zs space for plotting
zs_grid = floor(min(zs-thick_ice)):min(thick_ice)*dy:ceil(max(zs-thick_ice))+500;
zs_ages = zeros(length(zs_grid),length(x));
zs_ages(:,:) = NaN;
dUdx_ages = zeros(length(zs_grid),length(x));
dUdx_ages(:,:) = NaN;
dVdy_ages = zeros(length(zs_grid),length(x));
dVdy_ages(:,:) = NaN;

for n = 1:length(x)
        H_we_pt = interp1(x,H_we,x(n));
        zs_pt = interp1(dist,zs,x(n));
        zs_ages(:,n) = interp1(zs_pt-H_we_pt*rho_w/rho_i+zbar(2:end)*H_we_pt*rho_w/rho_i,ages(2:end,n),zs_grid);
        dUdx_ages(:,n) = interp1(zs_pt-H_we_pt*rho_w/rho_i+zbar(2:end)*H_we_pt*rho_w/rho_i,dUdx(2:end,n),zs_grid);
        dVdy_ages(:,n) = interp1(zs_pt-H_we_pt*rho_w/rho_i+zbar(2:end)*H_we_pt*rho_w/rho_i,dVdy(2:end,n),zs_grid);
end


%% Plot potential cores
potential_core_dists = [9e3,12.5e3,15e3,17e3];

fig1=figure('units','inches');
clear h
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 7 2.5],'color','white')

[ha,pos] = tight_subplot(1,4,0.01,0.17,0.1);
colors = [[0 0 0];[1 0 0];[0 0.5 0];[0 0 1]];

for i=1:length(potential_core_dists)
    axes(ha(i));  
    [~,ind1] = min(abs(potential_core_dists(i)-x));
    plot(ages(2:end,ind1)/1e3,(1-zbar(2:end))*H_ice(ind1),'color',colors(i,:),'linewidth',1.5);
    [~,ind2] = min(abs(potential_core_dists(i)-dist));
    hold on;
    [~,ind3] = min(sqrt((track1B.X1B_all-xtrack(ind2)).^2+(track1B.Y1B_all-ytrack(ind2)).^2));
    h(1) = plot(track1B.layages,track1B.layheights(:,ind3),'k.','markersize',15);
    interped = interp1((1-zbar(2:end))*H_ice(ind1),ages(2:end,ind1)/1e3,track1B.layheights(:,ind3));
    rmse=sqrt(nanmean((interped-track1B.layages).^2));
    set(gca,'ydir','reverse');
    set(gca,'fontsize',8,'fontname','Arial');
    if i==2
        xlabel('                                       Age (ka)','fontsize',8,'fontname','Arial');
    end
    set(gca,'ytick',[0:200:1200])
    if i==1
        ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    else
        set(gca,'yticklabel',[])
    end
    xlim([0,1200])
    ylim([0,1200])
    plot([0,1200],[track1B.thick(ind3),track1B.thick(ind3)],'k--','linewidth',1.5)
    set(gca,'ticklength',[0.035 0.07])
    text(200,100,sprintf('RMSE = %.1f ka',rmse),'fontsize',8,'fontname','Arial')
end

export_fig(fullfile(REPO_HOME,sprintf('figures/grinsted_model_cores_us%03d_bdot%02dmm.pdf',...
    us_constant*100,bdot_constant*1000)),...
    '-p0.01','-painters','-transparent');
close;

clear axes interped rmse ha pos fig1


%% Calculate RMSE between model and obs and save to file
rmse_dists = [5e3:dx:19e3];
sume = 0;
sumn = 0;
for i=1:length(rmse_dists)
    [~,ind1] = min(abs(rmse_dists(i)-x));
    [~,ind2] = min(abs(rmse_dists(i)-dist));
    [~,ind3] = min(sqrt((track1B.X1B_all-xtrack(ind2)).^2+(track1B.Y1B_all-ytrack(ind2)).^2));
    interped = interp1((1-zbar(2:end))*H_ice(ind1),ages(2:end,ind1)/1e3,track1B.layheights(:,ind3));
    nonnan = find(~(isnan(track1B.layages)) & ~(isnan(interped)));
    sume = sume+sum((interped(nonnan)-track1B.layages(nonnan)).^2);
    sumn = sumn+length(nonnan);
end
RMSE_TOTAL = sqrt(sume/sumn);
if save_rmse_file
    fprintf(fid,'%.2f %3f %.6f\n',us_constant,bdot_constant,RMSE_TOTAL);
end


%% Plot age structure from model and data for comparison
figure('units','inches');
pos = get(gcf,'pos');
cm = cubehelix(20,1.0,-1.5,1,1,[0.2,0.8],[0.2,1.0]);
[ha, pos] = tight_subplot(3,1,[0.00 0.00],[0.12 0.01],[0.15 0.01]);
set(gcf,'pos',[0 0 4 3],'color','white')
set(gca,'visible','off','box','off')

subplot(ha(1));
[ax h1 h2] = plotyy((x)/1e3, Us_back, (x)/1e3,bdot(:,1)*1000, 'plot');
hold(ax(1)); hold(ax(2));
set(h1,'LineStyle','--','Linewidth',1);
set(h2,'LineStyle','--','Linewidth',1,'Color','r');
set(ax(1),'ylim',[0 0.6],'ytick',[0:0.3:0.6],'fontsize',8,'fontname','Arial')
set(ax(2),'ylim',[-35 35],'ytick',[-30:30:30],'fontsize',8,'fontname','Arial')
set(ax(1),'ycolor','b','xlim',[(x(1))/1e3,(x(end))/1e3],'xtick',[-4:4:25]);
set(ax(2),'ycolor','r','xlim',[(x(1))/1e3,(x(end))/1e3],'xtick',[-4:4:25]);
plot(ax(1),(x_pts(1:end-1))/1e3,Us_pts(1:end-1),'ko','markerfacecolor','b','markersize',4);
ind = find(bdot_pts < 0);
plot(ax(2),(x_pts(ind))/1e3,bdot_pts(ind)*1000,'ko','markerfacecolor','r','markersize',4);
axes(ax(1)); ylabel('{\it{u_s}} (m/a)','fontsize',8,'fontname','Arial');
axes(ax(2)); ylabel({'{\it{b}}';'(mmWE/a)   '},'fontsize',8,'fontname','Arial');
text(22.2+dist_offset/1e3,11,'.','color','r','fontsize',8,'fontname','Arial');
pos = get(gca,'pos');
set(gca,'fontsize',8,'fontname','arial')
set(ax(1),'xticklabel',[])
set(ax(2),'xticklabel',[]);
set(gca,'ticklength',[0.02 0.05]);
text(-4+dist_offset/1e3,-10,'a','fontsize',10,'fontname','arial','fontweight','bold')

% Plot age structure from data
ax = subplot(ha(2));
set(gca,'fontsize',8,'fontname','Arial');
zdata = 800:1:2200;
if strcmp(flowline,'track4')
    distdata = track4.Distance_total4_all;
    filtdata = zeros(length(zdata),length(track4.Distance_total4_all));
    filtdata(:,:) = nan;
    ind = find(track4.z_out ~= 0);
    for i=1:length(track4.Distance_total4_all)
        [~,j] = min(abs(track4.Height4_all(i)-track4.thick(i)-zdata));
        filtdata(j:end,i) = interp1(track4.Height4_all(i)-track4.z_out(ind-1:end),track4.filtdata4_all(ind-1:end,i),zdata(j:end)); 
    end
elseif strcmp(flowline,'track1')
    distdata = track1B.Distance_total1B_all;
    filtdata = zeros(length(zdata),length(track1B.Distance_total1B_all));
    filtdata(:,:) = nan;
    ind = find(track1B.z_out ~= 0);
    for i=1:length(track1B.Distance_total1B_all)
        [~,j] = min(abs(track1B.Height1B_all(i)-track1B.thick(i)-zdata));
        filtdata(j:end,i) = interp1(track1B.Height1B_all(i)-track1B.z_out(ind-1:end),track1B.filtdata1B_all(ind-1:end,i),zdata(j:end)); 
    end
end
pcolor((fliplr(distdata)+dist_offset)/1e3,zdata,filtdata); hold on;
shading interp;
set(gca,'clim',[-0.0015 0.0015]);
colormap(flipud(colormap(bone)));
freezeColors;
if strcmp(flowline,'track4')
    for i = 1:length(track4.layages)
        if ~isnan(track4.layages(i))
            plot((x(end)-track4.Distance_total4_all+dist_offset)/1e3,track4.Height4_all-track4.layheights(i,:)',...
                'linewidth',3,'color',cm(floor(track4.layages(i)/50)+1,:));
        end
    end
    plot((track4.Distance_total4_all+dist_offset)/1e3,track4.Height4_all-track4.thick,'k','linewidth',2);   
elseif strcmp(flowline,'track1')
    for i = 1:length(track1B.layages)
        if ~isnan(track1B.layages(i))
            plot((fliplr(track1B.Distance_total1B_all)+dist_offset)/1e3,track1B.Height1B_all-track1B.layheights(i,:)',...
                'linewidth',3,'color',cm(floor(track1B.layages(i)/25)+1,:));
        else
            plot((fliplr(track1B.Distance_total1B_all)+dist_offset)/1e3,track1B.Height1B_all-track1B.layheights(i,:)',...
                'linewidth',1,'color','k');
        end
    end
    plot((fliplr(track1B.Distance_total1B_all)+dist_offset)/1e3,track1B.Height1B_all-track1B.thick,'k','linewidth',2);  
end
%potential_core_dists=[10e3];
%for i=1:length(potential_core_dists)
%    plot([potential_core_dists(i)/1e3,potential_core_dists(i)/1e3],[800,2200],'k--','linewidth',2);
%end
xlim([(x(1))/1e3,(x(end))/1e3])
ylim([800 2200])
set(gca,'ydir','normal');
ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
set(gca,'fontsize',8,'fontname','arial','ticklength',[0.02 0.05],'xtick',[-4:4:25])
pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2)+0.075 0.7 pos(4)+0.05])
text(-0.2+dist_offset/1e3,980,'1''','fontsize',8,'fontname','arial','fontweight','bold');
text(18.6+dist_offset/1e3,980,'1 ','fontsize',8,'fontname','arial','fontweight','bold');
text(-4+dist_offset/1e3,980,'b','fontweight','bold','fontsize',10,'fontname','arial');
set(gca,'xticklabel',[]);

ax = subplot(ha(3));
set(gca,'fontsize',8,'fontname','Arial');
pcolor((x)/1e3,zs_grid,zs_ages/1e3); hold on;
set(gca,'clim',[0 500]);
colormap(gca,cm);
xlabel('{\it{x}} (km)','fontsize',8,'fontname','Arial');

shading interp;
freezeColors;
plot((x)/1e3,zs_x,'k','linewidth',1.5);
hold on;
plot((x)/1e3,zs_x-H_ice,'k','linewidth',2);
%for i=1:length(potential_core_dists)
%    plot([potential_core_dists(i)/1e3,potential_core_dists(i)/1e3],[800,2200],'k:','linewidth',2);
%end
xlim([(x(1))/1e3,(x(end))/1e3]);
ylim([800 2200])
ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
set(gca,'fontsize',8,'fontname','arial','ticklength',[0.02 0.05],'xtick',[-4:4:25])
colormap(gca,cm); set(gca,'clim',[0,500]);
h = colorbar; set(h,'fontsize',8,'fontname','Arial','ytick',[0:100:500],'ticklength',[0.25 0.25]);
pos = get(h,'pos');
set(h,'pos',[pos(1) pos(2)+0.2 pos(3) pos(4)]);
set(get(h,'title'),'string','Age (ka)','fontsize',8,'fontname','arial');
pos = get(get(h,'title'),'position');
set(get(h,'title'),'position',[pos(1)+1 pos(2)+0.5 pos(3)])
pos2 = get(gca,'pos');
set(gca,'pos',[pos2(1) pos2(2)+0.005 0.7 pos2(4)+0.05-0.005]);
text(-4+dist_offset/1e3,980,'c','fontweight','bold','fontname','arial','fontsize',10);

subplot(ha(1));
pos = get(gca,'pos');
set(gca,'pos',[pos(1) pos(2)+0.15 0.7 pos(4)-0.15])

export_fig(fullfile(REPO_HOME,sprintf('figures/grinsted_model_us%03d_bdot%02dmm.pdf',...
    us_constant*100,bdot_constant*1000)),...
    '-p0.01','-zbuffer','-nofontswap');
close;

clear h cm ax pos colors


%% Plot upper packages
deviation_from_mean=1;

grinsted_layers_ages = [0e3:20e3:100e3];
% Determine ages of layers by comparing data to grinsted model
%ind = find((x >= 5e3) & (x <= 19e3));
%grinsted_layers_ages(2) = 1e3*nanmean(interp2(x(ind)-dist_offset,zs_grid,zs_ages(:,ind)/1e3,fliplr(track1B.Distance_total1B_all),track1B.layerpacks(1,:)));
%grinsted_layers_ages(3) = 1e3*nanmean(interp2(x(ind)-dist_offset,zs_grid,zs_ages(:,ind)/1e3,fliplr(track1B.Distance_total1B_all),track1B.layerpacks(2,:)));
%grinsted_layers_ages(4) = 1e3*nanmean(interp2(x(ind)-dist_offset,zs_grid,zs_ages(:,ind)/1e3,fliplr(track1B.Distance_total1B_all),track1B.layerpacks(3,:)));
%grinsted_layers_ages(5) = 1e3*nanmean(interp2(x(ind)-dist_offset,zs_grid,zs_ages(:,ind)/1e3,fliplr(track1B.Distance_total1B_all),track1B.layerpacks(4,:)));
%grinsted_layers_ages(6) = 1e3*nanmean(interp2(x(ind)-dist_offset,zs_grid,zs_ages(:,ind)/1e3,fliplr(track1B.Distance_total1B_all),track1B.layerpacks(5,:)));

grinsted_layers = zeros(length(x),length(grinsted_layers_ages));
grinsted_layers(:,1) = zs_x;
for j=2:length(grinsted_layers_ages)
    for i=1:length(x)
        nonnan = find(~(isnan(zs_ages(:,i))));
        [~,ind] = min(abs(zs_ages(:,i)-grinsted_layers_ages(j)));
        grinsted_layers(i,j) = zs_grid(ind);   
    end
    ind = find(abs(grinsted_layers(:,j)-grinsted_layers(:,1))<1);
    if length(ind) > 0
        grinsted_layers(ind(1):end,j) = NaN;
    end
end

grinsted_layers_diffs = zeros(length(x),length(grinsted_layers_ages)-1);
for i=1:length(grinsted_layers_ages)-1
    grinsted_layers_diffs(:,i) = grinsted_layers(:,i)-grinsted_layers(:,i+1);
    if deviation_from_mean
        grinsted_layers_diffs(:,i) = grinsted_layers_diffs(:,i) - nanmean(grinsted_layers_diffs(:,i));
    end
end

figure('units','inches');
coloptions = [[0,0,1];[0,1,1];[212/255,175/255,55/255];[1,0,0];[0,0,0]];

hold off;
pos = get(gcf,'pos');
set(gca,'XTickLabel','','YTickLabel','','box','off','visible','off');
set(gcf,'pos',[0 0 3.5 3],'color','white')
[ha, pos] = tight_subplot(2,1,[0.01 0.01],[0.25 0.02],[0.25 0.01]);

subplot(ha(2));
hold on;
box('on');
for i=1:size(grinsted_layers_diffs,2)
    plot((x)/1e3,grinsted_layers_diffs(:,i),'color',coloptions(i,:),'linewidth',3);
end
if strcmp(flowline,'track1')
    xlim([dist_offset/1e3,track1B.Distance_total1B_all(end)/1e3+dist_offset/1e3])
elseif strcmp(flowline,'track4')
    xlim([dist_offset/1e3,track4.Distance_total4_all(end)/1e3+dist_offset/1e3])
end
if deviation_from_mean
    set(ha(2),'ytick',[-100:50:100],'fontsize',8,'fontname','Arial')
    ylim([-120,120])
else
    ylim([0,300])
    set(ha(2),'ytick',[0:50:300],'fontsize',8,'fontname','Arial')    
end
ylabel({'Deviation from mean';'layer thickness (m)'},'fontsize',8,'fontname','Arial')
xlabel('{\it{x}} (km)','fontsize',8,'fontname','Arial')
set(gca,'ticklength',[0.02 0.05]);
% Legend
%rectangle('Position',[15.6,-10,3.8,120],'facecolor','w','edgecolor','k');
%plot([15.9,16.6],[96,96],'color',get(h(1),'color'),'linewidth',1.25)
%text(16.8,98,'0 - 1','fontsize',8,'fontname','arial');
%plot([15.9,16.6],[74,74],'color',get(h(2),'color'),'linewidth',1.25)
%text(16.8,76,'1 - 2','fontsize',8,'fontname','arial');
%plot([15.9,16.6],[52,52],'color',get(h(3),'color'),'linewidth',1.25)
%text(16.8,54,'2 - 3','fontsize',8,'fontname','arial');
%plot([15.9,16.6],[30,30],'color',get(h(4),'color'),'linewidth',1.25)
%text(16.8,32,'3 - 4','fontsize',8,'fontname','arial');
%plot([15.9,16.6],[8,8],'color',get(h(5),'color'),'linewidth',1.25)
%text(16.8,10,'4 - 5','fontsize',8,'fontname','arial');
text(0.5+dist_offset/1e3,100,'d','fontsize',9,'fontname','arial','fontweight','bold');

subplot(ha(1));
ylabel('Elevation (m)','fontsize',8,'fontname','Arial');
hold on;
box('on');
for i=1:size(grinsted_layers,2)-1
    nonnan = find((~(isnan(grinsted_layers(:,i)))) & (~(isnan(grinsted_layers(:,i+1)))));
    x_area1 = (x)/1e3;
    x_area = [x_area1(nonnan),fliplr(x_area1(nonnan))];
    y_area = [(grinsted_layers(nonnan,i))',(flipud(grinsted_layers(nonnan,i+1)))'];
    h(i) =  fill(x_area,y_area,coloptions(i,:),'facealpha',0.5,'edgecolor','none');
end
plot((x)/1e3,(zs_x-H_ice),'k','linewidth',1.5);
plot((x)/1e3,(zs_x),'k','linewidth',1.5);
plot((x)/1e3,grinsted_layers,'k','linewidth',1.5);
ylim([800,2200])
if strcmp(flowline,'track1')
    xlim([dist_offset/1e3,track1B.Distance_total1B_all(end)/1e3+dist_offset/1e3])
elseif strcmp(flowline,'track4')
    xlim([dist_offset/1e3,track4.Distance_total4_all(end)/1e3+dist_offset/1e3])
end
set(ha(1),'XTickLabel','');
set(ha(1),'layer','top','fontsize',8,'fontname','Arial');
rectangle('Position',[9.0+dist_offset/1e3,875,10.2,455],'facecolor','w','edgecolor','k');
rectangle('Position',[9.4+dist_offset/1e3,1185,0.6,80],'facecolor',coloptions(1,:),'edgecolor','k');
text(10.3+dist_offset/1e3,1240,sprintf('1 (%d ka)',round(grinsted_layers_ages(2)/1e3)),...
    'fontsize',8,'fontname','arial');
rectangle('Position',[9.4+dist_offset/1e3,1065,0.6,80],'facecolor',coloptions(2,:),'edgecolor','k');
text(10.3+dist_offset/1e3,1110,sprintf('2 (%d ka)',round(grinsted_layers_ages(3)/1e3)),...
    'fontsize',8,'fontname','arial');
rectangle('Position',[9.4+dist_offset/1e3,925,0.6,80],'facecolor',coloptions(3,:),'edgecolor','k');
text(10.3+dist_offset/1e3,980,sprintf('3 (%d ka)',round(grinsted_layers_ages(4)/1e3)),...
    'fontsize',8,'fontname','arial');
rectangle('Position',[14.2+dist_offset/1e3,1185,0.6,80],'facecolor',coloptions(4,:),'edgecolor','k');
text(15.1+dist_offset/1e3,1240,sprintf('4 (%d ka)',round(grinsted_layers_ages(5)/1e3)),...
    'fontsize',8,'fontname','arial');
rectangle('Position',[14.2+dist_offset/1e3,1065,0.6,80],'facecolor',coloptions(5,:),'edgecolor','k');
text(15.1+dist_offset/1e3,1110,sprintf('5 (%d ka)',round(grinsted_layers_ages(6)/1e3)),...
    'fontsize',8,'fontname','arial');
set(gca,'ticklength',[0.02 0.05]);
text(0.5+dist_offset/1e3,940,'b','fontweight','bold','fontsize',9,'fontname','arial');

clear h legend_h ha pos

export_fig(fullfile(REPO_HOME,sprintf('figures/grinsted_model_layers_us%03d_bdot%02dmm.pdf',...
    us_constant*100,bdot_constant*1000)),...
    '-opengl','-nofontswap','-r600');
close;


%% Plot particle trajectories
figure('units','inches');

hold off;
pos = get(gcf,'pos');
set(gca,'XTickLabel','','YTickLabel','','box','off');
set(gcf,'pos',[0 0 7 3],'color','white')

box('on');
plot(pt17(:,1)/1e3,pt17(:,3)); hold on;
plot(pt18(:,1)/1e3,pt18(:,3)); hold on;
plot(pt19(:,1)/1e3,pt19(:,3)); hold on;
plot(pt20(:,1)/1e3,pt20(:,3)); hold on;
plot(pt21(:,1)/1e3,pt21(:,3)); hold on;
plot(pt22(:,1)/1e3,pt22(:,3)); hold on;
plot(pt23(:,1)/1e3,pt23(:,3)); hold on;
plot(pt24(:,1)/1e3,pt24(:,3)); hold on;
plot(x/1e3,zs_x,'k','linewidth',1.5);
plot(x/1e3,zs_x-H_ice,'k','linewidth',1.5);
set(gca,'ydir','normal');
set(gca,'clim',[-8,8]);
xlabel('Distance (km)','fontsize',8,'fontname','Arial');
ylabel('Elevation (m)','fontsize',8,'fontname','Arial')
xlim([0,x(end)/1e3])

export_fig(fullfile(REPO_HOME,sprintf('figures/grinsted_model_points_us%03d_bdot%02dmm.pdf',...
    us_constant*100,bdot_constant*1000)),...
    '-p0.01','-zbuffer','-nofontswap');
close;

%% Save variables
if enhanced_bdot
    save(fullfile(REPO_HOME,sprintf('grinsted_matfiles/grinsted_us%03d_bdot%02dmm_enhancedbdot.mat',...
        us_constant*100,bdot_constant*1000)),...
    'RMSE_TOTAL','x','H_ice','H_we','grinsted_layers','grinsted_layers_ages','us_constant',...
    'bdot_constant','xtrack','ytrack','zs_grid','zs_ages','pt17','pt18','pt19','pt20','pt21','pt22',...
    'pt23','pt24','ages','Us','bdot','zbar','zs_x','dt','dx','dy','dist_offset','flowline');

else
    save(fullfile(REPO_HOME,sprintf('grinsted_matfiles/grinsted_us%03d_bdot%02dmm_highres.mat',...
    us_constant*100,bdot_constant*1000)),...
    'RMSE_TOTAL','x','H_ice','H_we','grinsted_layers','grinsted_layers_ages','us_constant',...
    'bdot_constant','xtrack','ytrack','zs_grid','zs_ages','pt17','pt18','pt19','pt20','pt21','pt22',...
    'pt23','pt24','ages','Us','bdot','zbar','zs_x','dt','dx','dy','dist_offset','flowline')
end

end
end

if save_rmse_file
    fclose(fid)
end

%% Plot RMSE
plot_rmse = 1;
if plot_rmse
    us_constants = [0:.1:1];
    bdot_constants = [0.005:0.001:0.014];
    rmse_grid = zeros([length(us_constants),length(bdot_constants)]);
    for i=1:length(bdot_constants)
        data = dlmread(fullfile(REPO_HOME,sprintf('RMSE/RMSE_bdot%02.0fmm.txt',bdot_constants(i)*1000)));
        rmse_grid(:,i) = data(:,3);
    end

    figure('units','inches');

    hold off;
    pos = get(gcf,'pos');
    set(gca,'XTickLabel','','YTickLabel','','box','off','fontsize',8,'fontname','Arial');
    set(gcf,'pos',[0 0 3.75 3],'color','white')


    imagesc(bdot_constants*1000,us_constants,rmse_grid);
    set(gca,'clim',[0,100])
    set(gca,'ydir','normal');
    cm = cubehelix(100,1.0,-1.5,1,1,[0.05,0.95],[0.0,1]);
    colormap((cm));
    h = colorbar; set(h,'fontsize',8,'fontname','Arial');
    set(get(h,'title'),'string','Misfit (ka)','fontsize',8,'fontname','Arial')
    set(h,'ticklength',[0.1 0.1])
    xlabel('Present-day accumulation rate ({\it{b_o}}; mmWE/a)','fontsize',8,'fontname','arial');
    set(gca,'xtick',bdot_constants*1000);
    set(gca,'ytick',us_constants);
    gpos = get(gca,'position');
    set(gca,'position',[gpos(1)+0.13 gpos(2) gpos(3)-0.01 gpos(4)])
    pos = get(h,'position');
    set(h,'position',[pos(1)+0.08 gpos(2) pos(3)*0.85 gpos(4)])
    ylabel({'Ratio of glacial to';' present-day velocities (\it{r})'},'fontsize',8,'fontname','arial');
    set(gca,'ticklength',[0.020,0.05])
    text(4.9,0.02,'a','fontweight','bold','fontsize',10,'fontname','arial');
    text(11.6,-0.13,'.','fontsize',8,'fontname','arial');

    export_fig(fullfile(REPO_HOME,'figures/rmse.pdf'),'-zbuffer','-nofontswap','-r600','-p0.01');
    close;
end
clear plot_rmse

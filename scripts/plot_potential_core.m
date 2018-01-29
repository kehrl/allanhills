% This script plots the depth-age scale for best-fit model combinations to
% determine how sensitive the age scale is to the chosen model. It also
% plots the depth of 1-Ma or 500-ka ice above the bed along the flowline.
%
% Laura Kehrl, UW, 10/01/2017

us030bdot10mm = load(fullfile(REPO_HOME,'grinsted_matfiles/grinsted_us030_bdot10mm_highres.mat'));
us050bdot09mm = load(fullfile(REPO_HOME,'grinsted_matfiles/grinsted_us050_bdot09mm_highres.mat'));
us060bdot08mm = load(fullfile(REPO_HOME,'grinsted_matfiles/grinsted_us060_bdot08mm_highres.mat'));
us080bdot07mm = load(fullfile(REPO_HOME,'grinsted_matfiles/grinsted_us080_bdot07mm_highres.mat'));

models = ['us080bdot07mm';'us060bdot08mm';'us050bdot09mm';'us030bdot10mm'];

%% Plot potential cores
%potential_core_dists = 10e3+dist_offset;
potential_core_dists = 16.0e3;

fig1=figure('units','inches');
clear h
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 3.6 2.92],'color','white')

colors = [[0 0 0];[1 0 0];[0 0 1];[0 0.5 0]];

for i=1:length(potential_core_dists)
    for j=1:size(models,1);
        disp(models(j,:));
        x = eval(sprintf('%s.x',models(j,:)));
        zbar = eval(sprintf('%s.zbar',models(j,:)));
        ages = eval(sprintf('%s.ages',models(j,:)));
        [~,ind1] = min(abs(potential_core_dists(i)-x));
        h(j)=plot(ages(2:end,ind1)/1e3,(1-zbar(2:end))*H_ice(ind1),'-','color',colors(j,:),'linewidth',1.25);
        [~,ind2] = min(abs(potential_core_dists(i)-dist));
        [~,ind3] = min(sqrt((track1B.X1B_all-xtrack(ind2)).^2+(track1B.Y1B_all-ytrack(ind2)).^2));
        hold on;
        interped = interp1((1-zbar(2:end))*H_ice(ind1),ages(2:end,ind1)/1e3,track1B.layheights(:,ind3));
        rmse=sqrt(nanmean((interped-track1B.layages).^2));
        age1 = 260e3;
        age2 = 240e3;
        [~,closeind] = min(abs(ages(:,ind1)-age1));
        if ages(closeind,ind1) < age1
            closeind1 = closeind-1;
            closeind2 = closeind;
        else
            closeind1 = closeind;
            closeind2 = closeind+1;
        end
        ma(i) = H_ice(ind1)-interp1(ages(closeind1:closeind2,ind1),(1-zbar(closeind1:closeind2))*H_ice(ind1),age1);
        [~,closeind] = min(abs(ages(:,ind1)-age2));
        if ages(closeind,ind1) < age2
            closeind1 = closeind-1;
            closeind2 = closeind;
        else
            closeind1 = closeind;
            closeind2 = closeind+1;
        end
        halfma(i) = H_ice(ind1)-interp1(ages(closeind1:closeind2,ind1),(1-zbar(closeind1:closeind2))*H_ice(ind1),age2);
        
        % Only works if it is monotonically increasing
        %ma = H_ice(ind1)-interp1(ages(3:end,ind1)/1e3,(1-zbar(2:end))*H_ice(ind1),751);
        %halfma = H_ice(ind1)-interp1(ages(2:end,ind1)/1e3,(1-zbar(2:end))*H_ice(ind1),749);
        rate(j) = 20/(halfma-ma);
    end
    plot(track1B.layages,track1B.layheights(:,ind3),'k.','markersize',15);    
    xlabel('Age (ka)','fontsize',8,'fontname','Arial');
    set(gca,'ytick',[0:200:1200])
    set(gca,'ydir','reverse');
    set(gca,'fontsize',8,'fontname','Arial');
    if i==1
        ylabel('Ice thickness (m)','fontsize',8,'fontname','Arial');
    else
        set(gca,'yticklabel',[])
    end
    xlim([0,1050])
    ylim([0,track1B.thick(ind3)])
    %plot([0,1200],[track1B.thick(ind3),track1B.thick(ind3)],'k--','linewidth',1.5)
    set(gca,'ticklength',[0.035 0.07])
end
%h=legend('{\it{b}}=7,   {\it{r}}=0.8','{\it{b}}=8,   {\it{r}}=0.6','{\it{b}}=9,   {\it{r}}=0.5','{\it{b}}=10, {\it{r}}=0.3','location','northwest');
rectangle('Position',[80,20,570,300],'facecolor','w','edgecolor','k');
plot([100,180],[65,65],'k','linewidth',1.25)
text(190,65,'{\it{b}} =   7 mmWE/a, {\it{r}} = 0.8','fontsize',8,'fontname','arial');
text(205,25,'.','fontsize',8,'fontname','arial');
plot([100,180],[135,135],'r','linewidth',1.25)
text(190,135,'{\it{b}} =   8 mmWE/a, {\it{r}} = 0.6','fontsize',8,'fontname','arial');
text(205,95,'.','fontsize',8,'fontname','arial');
plot([100,180],[205,205],'b','linewidth',1.25)
text(190,205,'{\it{b}} =   9 mmWE/a, {\it{r}} = 0.5','fontsize',8,'fontname','arial');
text(205,165,'.','fontsize',8,'fontname','arial');
plot([100,180],[275,275],'color',[0 0.5 0],'linewidth',1.25)
text(190,275,'{\it{b}} = 10 mmWE/a, {\it{r}} = 0.3','fontsize',8,'fontname','arial');
text(205,235,'.','fontsize',8,'fontname','arial');
text(25,1130,'b','fontsize',10,'fontname','arial','fontweight','bold');


fig2= figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 0.225 0.877],'color','white');
for i=1:length(potential_core_dists)
    for j=1:size(models,1);
        disp(models(j,:));
        x = eval(sprintf('%s.x',models(j,:)));
        zbar = eval(sprintf('%s.zbar',models(j,:)));
        ages = eval(sprintf('%s.ages',models(j,:)));
        [~,ind1] = min(abs(potential_core_dists(i)-x));
        plot(ages(2:end,ind1)/1e3,(1-zbar(2:end))*H_ice(ind1),'-','color',colors(j,:),'linewidth',1.25);
        [~,ind2] = min(abs(potential_core_dists(i)-dist));
        hold on;
        interped = interp1((1-zbar(2:end))*H_ice(ind1),ages(2:end,ind1)/1e3,track1B.layheights(:,ind3));
        rmse=sqrt(nanmean((interped-track1B.layages).^2));
    end
     [~,ind3] = min(sqrt((track1B.X1B_all-xtrack(ind2)).^2+(track1B.Y1B_all-ytrack(ind2)).^2));
    herrorbar(track1B.layages,track1B.layheights(:,ind3),ones([length(track1B.layages),1])*6,ones([length(track1B.layages),1])*6,'k.');    set(gca,'ytick',[0:200:1200])
    set(gca,'xtick',[100:100:200])
    set(gca,'ydir','reverse');
    set(gca,'fontsize',8,'fontname','Arial');
    xlim([60,225])
    ylim([450,1010])
    %plot([0,1200],[track1B.thick(ind3),track1B.thick(ind3)],'k--','linewidth',1.5)
    set(gca,'ticklength',[0.035 0.07])
end

inset(fig1,fig2,0.3,0.8);

export_fig(fullfile(REPO_HOME,'figures/potentialcore.pdf'),'-painters','-r600','-p0.01','-transparent','-nofontswap');
close;

%% Get depth of 1 Ma ice across basin for different model runs

age=1000e3;

fig = figure('units','inches');
pos = get(gcf,'pos');
set(gcf,'pos',[0 0 3.5 3],'color','white')
for j = 1:size(models,1)
    x = eval(sprintf('%s.x',models(j,:)));
    zbar = eval(sprintf('%s.zbar',models(j,:)));
    ages = eval(sprintf('%s.ages',models(j,:)));
    for i=1:length(x)
        if max(ages(:,i)) < age
            z1ma(i) = NaN;
        else
            [junk,ind1] = min(abs((ages(:,i) - age)));
            if ages(ind1,i) < age
                ind2 = ind1-1;
            else
                ind2 = ind1;
                ind1 = ind2+1;
            end
            if ages(ind2,i)==ages(ind1,i)
                z1ma(i) = (zbar(ind2)+zbar(ind1))/2*H_ice(i);
            else
                z1ma(i) = interp1(ages(ind2:ind1,i),(zbar(ind2:ind1))*H_ice(i),age);
            end
        end
    end
    plot((x)/1e3,z1ma,'-','color',colors(j,:),'linewidth',1.25);
    hold on;
end
xlabel('{\it{x}} (km)','fontsize',8,'fontname','arial')
set(gca,'xtick',[-4:4:20]);
xlim([0,20])
legend('{\it{b}}=7,   {\it{r}}=0.8','{\it{b}}=8,   {\it{r}}=0.6','{\it{b}}=9,   {\it{r}}=0.5','{\it{b}}=10, {\it{r}}=0.3','location','southeast');
set(gca,'fontsize',8,'fontname','arial')
if age == 1e6
    ylabel('Height of 1-Ma ice above bed (m)','fontsize',8,'fontname','arial');
    text(-4+dist_offset/1e3,33,'b','fontsize',8,'fontname','arial','fontweight','bold');
    export_fig('figures/depthto1ma.pdf','-painters','-r600','-p0.01','-transparent','-nofontswap');
else
    ylim([0,75])
    ylabel('Height of 500-ka ice above bed (m)','fontsize',8,'fontname','arial');
    text(-4+dist_offset/1e3,71,'a','fontsize',8,'fontname','arial','fontweight','bold');
    export_fig(fullfile(REPO_HOME,'figures/depthto500ka.pdf'),'-painters','-r600','-p0.01','-transparent','-nofontswap');
end

close;
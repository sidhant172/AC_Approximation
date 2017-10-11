%% Worst-case error in active power flow on each branch
% Based on optimal linearizations code (pscc2018):

clear
clc
close all


% load case14_errorbound_results
load case24_errorbound_results
% load case30_errorbound_results

var_pct = 0.1:0.1:0.4;

Sbase = 100;

ax = figure;

fsize = 14;

var_pct_plot = 0.4;
var_pct_idx = find(var_pct == var_pct_plot,1);

nanU = isnan(Pflow_diff_ub_full(:,var_pct_idx));
nanL = isnan(Pflow_diff_lb_full(:,var_pct_idx));
Pflow_diff_ub_full(nanU,var_pct_idx) = 0;
Pflow_diff_lb_full(nanL,var_pct_idx) = 0;

range = Pflow_diff_ub_full(:,var_pct_idx) - Pflow_diff_lb_full(:,var_pct_idx);
weight_diff = [Pflow_diff_lb_full(:,var_pct_idx) range];
bh = bar((1:size(mpc.branch,1)).',Sbase*weight_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none');%,'barwidth',0.4);
set(bh(2),'FaceColor',[255 204 204]/255,'EdgeColor','none');%,'barwidth',0.4);
hold on
for i=1:size(mpc.branch,1)
    plot([i-0.4 i+0.4],Pflow_diff_ub_full(i,var_pct_idx)*ones(2,1)*Sbase,'r-','LineWidth',2);
    plot([i-0.4 i+0.4],Pflow_diff_lb_full(i,var_pct_idx)*ones(2,1)*Sbase,'r--','LineWidth',1);
end


set(gca, 'YGrid', 'on', 'XGrid', 'off', 'FontSize',fsize,'XTick',1:size(mpc.branch,1))
xlabel('Line Index','FontSize',fsize);
ylabel('Possible Range of Worst-Case Power Flow Error (MW)','FontSize',fsize);
% title(['Power Flow Error at Each Bus: IEEE ' num2str(nbus) '-Bus System, ' num2str(var_pct(var_pct_idx)*100) '% Variation'],'FontSize',fsize);
title([num2str(var_pct(var_pct_idx)*100) '% Allowed Variation in Power Injections'],'FontSize',fsize);
a = axis;
axis([0.5 size(mpc.branch,1)+0.5 0 a(4)]);

set(gcf,'PaperPosition',[0.2500    2.5000    8.0000    6.0000]);

ax = figure;
var_pct_idx = find(var_pct == var_pct_plot,1);
range = Pflow_diff_ub_full(:,var_pct_idx) - Pflow_diff_lb_full(:,var_pct_idx);
weight_diff = [Pflow_diff_lb_full(:,var_pct_idx) range];
bh = bar((1:size(mpc.branch,1)).',Sbase*weight_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none');%,'barwidth',0.4);
set(bh(2),'FaceColor',[255 204 204]/255,'EdgeColor','none');%,'barwidth',0.4);
hold on
for i=1:size(mpc.branch,1)
    plot([i-0.4 i+0.4],Pflow_diff_ub_full(i,var_pct_idx)*ones(2,1)*Sbase,'r-','LineWidth',2);
    plot([i-0.4 i+0.4],Pflow_diff_lb_full(i,var_pct_idx)*ones(2,1)*Sbase,'r--','LineWidth',1);
end

set(gca, 'YGrid', 'on', 'XGrid', 'off', 'FontSize',fsize,'XTick',1:size(mpc.branch,1))
xlabel('Line Index','FontSize',fsize);
ylabel('Possible Range of Worst-Case Power Flow Error (MW)','FontSize',fsize);
% title(['Power Injection Error at Each Bus: IEEE ' num2str(nbus) '-Bus System, ' num2str(var_pct(var_pct_idx)*100) '% Variation'],'FontSize',fsize);
title([num2str(var_pct(var_pct_idx)*100) '% Allowed Variation in Power Injections'],'FontSize',fsize);
axis([0.5 size(mpc.branch,1)+0.5 0 a(4)]);

set(gcf,'PaperPosition',[0.2500    2.5000    8.0000    6.0000]);

% title(['Power Injection Error at Each Bus: IEEE ' num2str(nbus) '-Bus System'],'FontSize',fsize);
% a = axis;
% axis([a(1) a(2) 0 a(4)])


%% Worst-case error in reactive power flow on each branch
% Plot error (MVAr) at all branches for a given allowed variation
% Show gap between upper and lower bounds (fill in with light red)
%
ax = figure;

nanU = isnan(Qflow_diff_ub_full(:,var_pct_idx));
nanL = isnan(Qflow_diff_lb_full(:,var_pct_idx));
Qflow_diff_ub_full(nanU,var_pct_idx) = 0;
Qflow_diff_lb_full(nanL,var_pct_idx) = 0;

range = Qflow_diff_ub_full(:,var_pct_idx) - Qflow_diff_lb_full(:,var_pct_idx);
weight_diff = [Qflow_diff_lb_full(:,var_pct_idx) range];
bh = bar((1:size(mpc.branch,1)).',Sbase*weight_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none');%,'barwidth',0.4);
set(bh(2),'FaceColor',[255 204 204]/255,'EdgeColor','none');%,'barwidth',0.4);
hold on
for i=1:size(mpc.branch,1)
    plot([i-0.4 i+0.4],Qflow_diff_ub_full(i,var_pct_idx)*ones(2,1)*Sbase,'r-','LineWidth',2);
    plot([i-0.4 i+0.4],Qflow_diff_lb_full(i,var_pct_idx)*ones(2,1)*Sbase,'r--','LineWidth',1);
end


set(gca, 'YGrid', 'on', 'XGrid', 'off', 'FontSize',fsize,'XTick',1:size(mpc.branch,1))
xlabel('Line Index','FontSize',fsize);
ylabel('Possible Range of Worst-Case Power Flow Error (MVAr)','FontSize',fsize);
% title(['Power Flow Error at Each Bus: IEEE ' num2str(nbus) '-Bus System, ' num2str(var_pct(var_pct_idx)*100) '% Variation'],'FontSize',fsize);
title([num2str(var_pct(var_pct_idx)*100) '% Allowed Variation in Power Injections'],'FontSize',fsize);
a = axis;
axis([0.5 size(mpc.branch,1)+0.5 0 a(4)]);

set(gcf,'PaperPosition',[0.2500    2.5000    8.0000    6.0000]);

ax = figure;
var_pct_idx = find(var_pct == var_pct_plot,1);
range = Qflow_diff_ub_full(:,var_pct_idx) - Qflow_diff_lb_full(:,var_pct_idx);
weight_diff = [Qflow_diff_lb_full(:,var_pct_idx) range];
bh = bar((1:size(mpc.branch,1)).',Sbase*weight_diff,'stacked');
set(bh(1),'FaceColor','none','EdgeColor','none');%,'barwidth',0.4);
set(bh(2),'FaceColor',[255 204 204]/255,'EdgeColor','none');%,'barwidth',0.4);
hold on
for i=1:size(mpc.branch,1)
    plot([i-0.4 i+0.4],Qflow_diff_ub_full(i,var_pct_idx)*ones(2,1)*Sbase,'r-','LineWidth',2);
    plot([i-0.4 i+0.4],Qflow_diff_lb_full(i,var_pct_idx)*ones(2,1)*Sbase,'r--','LineWidth',1);
end

set(gca, 'YGrid', 'on', 'XGrid', 'off', 'FontSize',fsize,'XTick',1:size(mpc.branch,1))
xlabel('Line Index','FontSize',fsize);
ylabel('Possible Range of Worst-Case Power Flow Error (MVAr)','FontSize',fsize);
% title(['Power Injection Error at Each Bus: IEEE ' num2str(nbus) '-Bus System, ' num2str(var_pct(var_pct_idx)*100) '% Variation'],'FontSize',fsize);
title([num2str(var_pct(var_pct_idx)*100) '% Allowed Variation in Power Injections'],'FontSize',fsize);
axis([0.5 size(mpc.branch,1)+0.5 0 a(4)]);

set(gcf,'PaperPosition',[0.2500    2.5000    8.0000    6.0000]);

    

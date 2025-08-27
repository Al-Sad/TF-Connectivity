% Copyright (c) 2025 Mohammad Al-Sa'd
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%
% Email: ext-mohammad.al-sad@hus.fi, mohammad.al-sad@helsinki.fi
%        mohammad.al-sad@tuni.fi, alsad.mohamed@gmail.com
%
% The following reference should be cited whenever this script is used:
% Al-Sa'd, M., & Boashash, B. (2025). Time-Frequency Functional Brain 
% Networks: Application to Newborn Seizure Analysis. 
% IntechOpen. doi: 10.5772/intechopen.1011395
%
% Last Modification: 27-August-2025
%
% Description:
% This demo script generates figures 2 and 3 that illusrate the group 
% comparisons using the traditional and proposed time-frequency FBN measures.

%% Initialization
clear; close all; clc;
addpath(genpath('Functions'));

%% Parameters
opt_str = 'Conventional'; % 'Conventional' or 'Time Frequency'

%% Brain network processing
dataset_path = fullfile("Brain Networks",opt_str);
M = length(dir(dataset_path + "/*.mat"));
background_networks = [];
seizure_networks    = [];
healthy_networks    = [];
for m = 1:M  % iterate through the patients
    %%% Load the brain networks and check the subject status (epileptic or healthy)
    fprintf('Processing the brain network of subject %d/%d\n',m,M);
    temp = load(fullfile(dataset_path,"network" + m + ".mat"));
    subject_network = cat(5,temp.pear_network,temp.spear_network, ...
        temp.plv_network,temp.iplv_network,temp.ciplv_network, ...
        temp.pli_network,temp.wpli_network,temp.dwpli_network);
    subject_status  = sum(temp.eeg_mask) > 0;
    %%% Extract seizure and non-seizure networks
    seizure     = subject_network(temp.eeg_mask,:,:,:,:);
    non_seizure = subject_network(~temp.eeg_mask,:,:,:,:);
    clear subject_network;
    %%% Remove temporal outliers
    for i = 1:size(seizure,2) % frequency bands
        for p = 1:size(seizure,3) % electrode p
            for q = 1:size(seizure,4) % electrode q
                for j = 1:size(seizure,5) % connectivity measures
                    temp1 = squeeze(seizure(:,i,p,q,j));
                    temp2 = squeeze(non_seizure(:,i,p,q,j));
                    [~,I] = rmoutliers(temp1(:));
                    temp1(I) = nan;
                    [~,I] = rmoutliers(temp2(:));
                    temp2(I) = nan;
                    seizure(:,i,p,q,j) = temp1;
                    non_seizure(:,i,p,q,j) = temp2;
                end
            end
        end
    end
    %%% Average the seizure and non-seizure networks across time
    seizure     = squeeze(mean(seizure,'omitnan'));
    non_seizure = squeeze(mean(non_seizure,'omitnan'));
    %%% Extract seizure, background, and healthy networks
    if subject_status
        background_networks = cat(5,background_networks,non_seizure);
        seizure_networks    = cat(5,seizure_networks,seizure);
    else
        healthy_networks = cat(5,healthy_networks,non_seizure);
    end
end

%% Average the brain networks of each group
healthy_avg_network = squeeze(mean(healthy_networks,5,'omitnan'));
background_avg_network = squeeze(mean(background_networks,5,'omitnan'));
seizure_avg_network = squeeze(mean(seizure_networks,5,'omitnan'));
I = permute(logical(repmat(eye(size(healthy_avg_network,2)), ...
    [1,1,size(healthy_avg_network,[1 4])])),[3 1 2 4]);
healthy_avg_network(I) = nan;
background_avg_network(I) = nan;
seizure_avg_network(I) = nan;

%% Multiple comparison test
multi_compare_test = cell(size(seizure_avg_network,1),size(seizure_avg_network,4));
for k = 1:size(seizure_avg_network,1)
    for j = 1:size(seizure_avg_network,4)
        x1 = squeeze(healthy_avg_network(k,:,:,j));
        x2 = squeeze(background_avg_network(k,:,:,j));
        x3 = squeeze(seizure_avg_network(k,:,:,j));
        x1(isnan(x1)) = [];
        x2(isnan(x2)) = [];
        x3(isnan(x3)) = [];
        [p,t,stats] = kruskalwallis([x1, x2, x3], ...
            [ones(size(x1)), 2*ones(size(x2)), 3*ones(size(x3))],'off');
        multi_compare_test{k,j} = multcompare(stats,"CriticalValueType", ...
            "bonferroni","Display","off");
    end
end

%% Plotting
switch opt_str
    case 'Conventional'
        grpLabels = {'Pearson','Spearman','PLV','iPLV','ciPLV','PLI','wPLI','dwPLI'};
    case 'Time Frequency'
        grpLabels = {'TF-Pearson','TF-Spearman','TF-PLV',...
            'TF-iPLV','TF-ciPLV','TF-PLI','TF-wPLI','TF-dwPLI'};
end
Color = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.8500 0.3250 0.0980]};
Labels = {'$\delta$','$\theta$','$\alpha$','$\beta$'};

figure('Color','w','Position',[50 -500 1200 1200]);
tiledlayout(size(seizure_avg_network,1),length(grpLabels), ...
    'TileSpacing','tight','Padding','compact');
for k = 1:size(seizure_avg_network,1)
    for j = 1:size(seizure_avg_network,4)
        t = nexttile; hold(t);
        x1 = squeeze(healthy_avg_network(k,:,:,j));
        x2 = squeeze(background_avg_network(k,:,:,j));
        x3 = squeeze(seizure_avg_network(k,:,:,j));
        x1(isnan(x1)) = [];
        x2(isnan(x2)) = [];
        x3(isnan(x3)) = [];
        data = {x1, x2, x3};
        violinPlot(data,'color',Color,'distWidth',1,'histOpt',1,'showMM',6)
        Max1 = min(max([x1,x2,x3])+max([std(x1) std(x2) std(x3)]),1);
        Min1 = max(min([x1,x2,x3])-max([std(x1) std(x2) std(x3)]),0);
        Max2 = min(max([x1,x2,x3])+0.50.*max([std(x1) std(x2) std(x3)]),1);
        Min2 = max(min([x1,x2,x3])-0.50.*max([std(x1) std(x2) std(x3)]),0);
        Max3 = min(max([x1,x2,x3])+0.75.*max([std(x1) std(x2) std(x3)]),1);
        Min3 = max(min([x1,x2,x3])-0.75.*max([std(x1) std(x2) std(x3)]),0);
        if multi_compare_test{k,j}(1,end) < 0.001
            plot([1 2], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(1.3, Min3, '*k','MarkerSize',6);
            plot(1.5, Min3, '*k','MarkerSize',6);
            plot(1.7, Min3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(1,end) < 0.01
            plot([1 2], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(1.4, Min3, '*k','MarkerSize',6);
            plot(1.6, Min3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(1,end) < 0.05
            plot([1 2], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(1.5, Min3, '*k','MarkerSize',6);
        end
        if multi_compare_test{k,j}(3,end) < 0.001
            plot([2 3], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(2.3, Min3, '*k','MarkerSize',6);
            plot(2.5, Min3, '*k','MarkerSize',6);
            plot(2.7, Min3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(3,end) < 0.01
            plot([2 3], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(2.4, Min3, '*k','MarkerSize',6);
            plot(2.6, Min3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(3,end) < 0.05
            plot([2 3], [Min2 Min2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(2.5, Min3, '*k','MarkerSize',6);
        end
        if multi_compare_test{k,j}(2,end) < 0.001
            plot([1 3], [Max2 Max2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(1.8, Max3, '*k','MarkerSize',6);
            plot(2, Max3, '*k','MarkerSize',6);
            plot(2.2, Max3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(2,end) < 0.01
            plot([1 3], [Max2 Max2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(1.9, Max3, '*k','MarkerSize',6);
            plot(2.1, Max3, '*k','MarkerSize',6);
        elseif multi_compare_test{k,j}(2,end) < 0.05
            plot([1 3], [Max2 Max2], '-k', 'LineWidth', 1.5, 'Marker','|');
            plot(2, Max3, '*k','MarkerSize',6);
        end
        box on; grid on; axis([0.4 3.6 Min1 Max1]);
        if k == size(seizure_avg_network,1)
            set(gca,'XTick',1:3,'XTickLabel',{'1','2','3'}, ...
                'TickLabelInterpreter','latex','FontSize',12);
            xlabel(grpLabels{j},'Interpreter','latex','FontSize',12);            
        else
            set(gca,'XTick',1:3,'XTickLabel','','TickLabelInterpreter','latex','FontSize',12);
        end
        if j == 1
            ylabel([Labels{k} '-Band'],'Interpreter','latex','fontsize',20);
            if k == 1
                legend({'\,\,Healthy','\,\,Non-Seizure','\,\,Seizure'},'FontSize',11, ...
                    'Interpreter','latex','EdgeColor','none','Color','none', ...
                    'Orientation','horizontal','Position',[0.06252 0.963 0.28374 0.021327]);
                axis([0.4 3.6 Min1-0.1 Max1]);
            end
        end

    end
end
set(gcf,'PaperUnits','inches','Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    if ~isfolder("Results"), mkdir("Results"); end
    print(1,fullfile("Results","comparison_" + lower(opt_str)),'-dpdf','-r400');
end

%
% FILE NAME: gen_TNBCprj_tINITmodels.m
%
%
%clear;

%% Clone and add needed packages and libraries to Matlab path
% cd /Users/fariba/Documents/Metastatic-TNBC/github
% system('git clone https://github.com/SysBioChalmers/RAVEN.git');
% cd RAVEN
% system('git checkout cc3f6e3')
% cd ..
% system('git clone https://github.com/SysBioChalmers/Human-GEM.git');
% cd Human-GEM
% system('git checkout 74de28c')
%cd ..


%Add them to path
addpath(genpath('/Users/fariba/Documents/GitHub/RAVEN'))
addpath(genpath('/Users/fariba/Documents/GitHub/human-GEM'))
addpath(genpath('/Users/fariba/Documents/Metastatic-TNBC/scripts'))


%% Load necessary files and data
cd /Users/fariba/Documents/Metastatic-TNBC/data
load('HumanGEM.mat');
length(ihuman.rxns)
% load metabolic tasks file
%taskStruct = parseTaskList('metabolicTasks_Essential.xlsx');
%load arrayData
load('bC2_refbGTex_nonlog_median_ArrayData.mat');

% Run some preliminary steps that will allow us to skip some pre-processing
% steps in the tINIT algorithm, greatly reducing the overall run time.
% refModel = ihuman;
% [~, deletedDeadEndRxns] = simplifyModel(refModel,true,false,true,true,true);
% cModel = removeReactions(refModel,deletedDeadEndRxns,false,true);
% [taskReport, essentialRxnMat] = checkTasks(cModel,[],true,true,true,taskStruct);
% essentialRxnsForTasks = cModel.rxns(any(essentialRxnMat,2));
%save these inputs for later use
%save('HumanGEMv1_3_2_fastInit_prepFiles','cModel','deletedDeadEndRxns','essentialRxnMat','taskReport','taskStruct')
load ('HumanGEMv1_3_2_fastInit_prepFiles')

%% Run tINIT algorithm
% the following script will be divided to several scripts to be run on the
% cluster
for i = 1:length(arrayData)
clear init_model;
    % add these shortcut fields to arrayData structure
    arrayData{i}.deletedDeadEndRxns = deletedDeadEndRxns;
    arrayData{i}.taskReport = taskReport;
    arrayData{i}.essentialRxnMat = essentialRxnMat;
    % run getINITModel function (use the special "fast" function version)
    [init_model, metProduction, essentialRxnsForTasks, addedRxnsForTasks, deletedDeadEndRxns, deletedRxnsInINIT, taskReport] = ...
        getINITModel2_fast(ihuman,char(arrayData{i}.tissues),[],[],arrayData{i},[],true,[],true,false,taskStruct);
    %save tINIT model
    cd('Models/tINIT_models_Median')
save(['tINIT_model_' num2str(i) '.mat'],'init_model')
cd('../..')
end 


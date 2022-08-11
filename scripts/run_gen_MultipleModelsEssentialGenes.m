%
% FILE NAME:    run_gen_MultipleModelsEssentialGenes.m
%
% DATE CREATED: 2020-01-27
%
% PROGRAMMER:   Fariba Roshanzamir
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
%  
% PURPOSE:  To automate doing  gene essentiality analysis of reconstructed
%           models in TNBCmetastatic project which includes Primary Tumor,
%           Adjacent Normal , Healthy and also TNBC metstatic models for 
%           different tissues.
%  
%
% REQUIRED VARIABLES:
%
%   inputModelsFile     The address to the structure containing simplified_tINIT models
%                       (or information necessary to regenerate tINIT models)
%                       for which gene essentiality will be evaluated.  
%
%   refModelFile        reference model from which the tINIT models were
%                       generated.
%
%   taskFile            metabolic task structure
%
% Fariba Roshanzamir 2020-01-27


%%
% Script to identify genes essential for metabolic tasks
% NOTE: these functions require that you have a working installation of
% gurobi and also raven and cobratoolbox
%  
% cd /Library/gurobi811/mac64/matlab
% gurobi_setup
% savepath
setRavenSolver('gurobi');


%clear
%___________  add needed toolboxes and functions to the working path  _____
%addpath /Applications/CPLEX_Studio128/cplex/matlab;
% addpath /Users/fariba/Documents/MATLAB/Human-GEM-devel20190722;
% addpath /Users/fariba/Documents/MATLAB/RAVEN-compareMultipleModels20190722;
% addpath /Users/fariba/Documents/MATLAB/cobratoolbox;
% addpath /Users/fariba/Documents/BrCprojectNewdataset/matlabScripts;
% addpath /Library/gurobi811/mac64/matlab

addpath(genpath('/Users/fariba/Documents/GitHub/RAVEN'))
addpath(genpath('/Users/fariba/Documents/GitHub/human-GEM'))
addpath(genpath('/Users/fariba/Documents/GitHub/cobratoolbox'))
addpath /Library/gurobi811/mac64/matlab
addpath '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Scripts_Matlab'

%________________  prepare REQUIRED VARIABLES  ____________________________
% Put all needed files in one folder and address to them
filepath = '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Analysis/';
x = dir(filepath);
Modelfiles = ({x.name})';
Modelfiles(cellfun(@isempty, regexp(Modelfiles,'TNBCmetastaticProj_INITmodels.mat$'))) = [];
modelsfileAddress = strcat(filepath, Modelfiles)

taskfile = ({x.name})';
taskfile(cellfun(@isempty, regexp(taskfile,'metabolicTasks_Essential.xlsx$'))) = [];
taskfileAddress =  strcat(filepath, taskfile)

% load the reference humanGEM model
load('humanGEM.mat'); % loads humanGEM model as variable "ihuman"

%________________  gene essentiality analysis  ____________________________
% run gen_MultipleModelsEssentialGenes function
inputModelsFile = modelsfileAddress;
refmodel = ihuman;
taskfile = taskfileAddress;
essGenes = gen_MultipleModelsEssentialGenes(inputModelsFile,refmodel,taskfile);


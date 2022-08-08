%
% FILE NAME:    gen_MultipleModelsEssentialGenes.m
% 
% DATE CREATED: 2020-01-27
% 
% PROGRAMMER:   Fariba Roshanzamir
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: To generate results of gene essentiality analysis for all
%          reconstructed models in metastatic TNBC breast cancer project
%
% Input:
%
%   inputModelsFile    The address to the structure containing simplified_tINIT models 
%                      for which gene essentiality will be evaluated.
%                      
%
%   refModel           reference model from which the tINIT models were
%                      generated.
%
%   taskfile           metabolic task file address
%
% Output:
%
%   essGenes          results structure with the following fields:
%       taskList      list of metabolic tasks that were tested
%       tissues       list of tissues (model IDs) corresponding to each model
%       geneList      cell array of the list of genes from each model
%       essentialGenes   cell array with one entry per model, where each
%                        entry is a logical matrix with rows corresponding
%                        to genes (in geneList) and columns to tasks (in
%                        taskList). Entries in the matrix are true when a
%                        gene is essential for a task, and false otherwise.

%
% Fariba Roshanzamir 2020-01-27


function essGenes = gen_MultipleModelsEssentialGenes(inputModelsFile,refmodel,taskfile)

% load the reference humanGEM model
% loads humanGEM model as variable "ihuman"
% assign refModel
refModel = refmodel;
% convert refModel genes from Ensembl IDs to gene abbreviations
[grRules,genes,rxnGeneMat] = translateGrRules(refModel.grRules,'Name','ENSG');
refModel.grRules = grRules;
refModel.genes = genes;
refModel.rxnGeneMat = rxnGeneMat;

% load the tINIT models
% these models are simplified so we need to bring back Boundary Metabolites
% so I reproduce them here to be used as input for getTaskEssentialGenes
% which will use "checkTasksGenes" function 
simpINITmodels = load(inputModelsFile{1});
simpINITmodels = simpINITmodels.TNBCmetastaticProj_INITmodels; %this isn't a good way for automate programming

for i=1:length(simpINITmodels)
%Add Boundary Metabolites Because Exchange metabolites should normally not be removed from the model when using checkTasks
comparingINIT_Models.model{i,1} = addBoundaryMets(simpINITmodels{1,i});
comparingINIT_Models.id{i,1} = simpINITmodels{1,i}.id;
end

% load the metabolic tasks
taskStruct = parseTaskList(taskfile{1});
% Identify genes essential for different tasks in different tINIT models.
%taskStruct =taskStruct(57); %  only check biomass production task for now (last task in the list)
taskStruct =taskStruct; %  only check biomass production task for now (last task in the list)

% get essential genes for each model and task
essGenes = getTaskEssentialGenes(comparingINIT_Models, refModel, taskStruct);

modelfile =inputModelsFile{1};

%filename = regexprep(modelfile, '\.mat', '_essGenes.mat');
filename = regexprep(modelfile, '\.mat', '_essGenes_AllTasks.mat');

save(filename, 'essGenes'); 
end

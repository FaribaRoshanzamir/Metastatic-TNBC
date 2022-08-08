%
% FILE NAME:            compareGeneTaskEssentiality.m
% CHANGED FILE NAME:    compareGeneTaskEssentiality_Results.m
% 
% DATE CREATED: 2019-02-14
% DATE EDITED:  2020-01-28 Fariba Roshanzamir
% 
% NOTE: *** This script has been Edited to be used in TNBC project. ***
%
% PROGRAMMER:   Jonathan Robinson
%               Department of Biology and Biological Engineering
%               Chalmers University of Technology
% 
% PURPOSE: This script combines and visualizes differences between gene
%          essentiality results for different models. Gene essentiality
%          results must first be obtained using the "checkTasksGenes"
%          function (with the "getEssential" option set to TRUE).
%
%
% REQUIRED VARIABLES:
%
%   models               A cell array of model structures, where each model
%                        has a different "id" field that will be used to
%                        differentiate the models.
%
%   essentialGeneArray   A cell array of essential gene matrices
%                        corresponding to each model. Each essential gene
%                        matrix is an MxN matrix with the essential genes
%                        (M) for each task (N). An element is true if the
%                        corresponding gene is essential in the
%                        corresponding task.
%
%                        Note: the different essential gene matrices for
%                        the different models do not need to have the same
%                        number of rows (genes), since they will be aligned
%                        automatically. However, they do need to have the
%                        same number and ordering of columns (metabolic
%                        tasks).
%
%   taskNames            A list of names or IDs of each task corresponding
%                        to the columns of the essentialGene matrix.
%
%
% Jonathan Robinson, 2019-02-14
% 

filepath = '/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Analysis/';
x = dir(filepath);
Modelfiles = ({x.name})';
Modelfiles(cellfun(@isempty, regexp(Modelfiles,'TNBCmetastaticProj_INITmodels.mat$'))) = [];
modelsfileAddress = strcat(filepath, Modelfiles)
%
load(modelsfileAddress{1})
models = TNBCmetastaticProj_INITmodels;
modelIDs = cellfun(@(m) {m.id},models);


if numel(unique(modelIDs)) < numel(models)
    error('All model IDs must be different.');
end


Modelfiles = ({x.name})';
Modelfiles(cellfun(@isempty, regexp(Modelfiles,'TNBCmetastaticProj_INITmodels_essGenes_AllTasks.mat$'))) = [];
essGenesfileAddress = strcat(filepath, Modelfiles)
eGenes = load(essGenesfileAddress{1})
essentialGeneArray = eGenes.essGenes;
taskNames = eGenes.essGenes.taskList;
genes = {};


for i = 1:numel(models)
    [grRules,genes,rxnGeneMat] = translateGrRules(models{i}.grRules,'Name','ENSG');
    models{i}.grRules = grRules;
    models{i}.genes = genes;
    models{i}.rxnGeneMat = rxnGeneMat;
end

for i = 1:numel(models)
    if numel(models{i}.genes) ~= size(essentialGeneArray.essentialGenes{i},1)
        error('The rows of each matrix in essentialGeneArray must equal the number of genes in the corresponding model.');
    end
    genes = union(genes,models{i}.genes);
end
length(genes)
eMat = zeros(numel(genes),numel(taskNames),numel(models));
for i = 1:numel(models)
    [isPresent,geneInd] = ismember(genes,models{i}.genes);
    eMat(isPresent,:,i) = essentialGeneArray.essentialGenes{i}(geneInd(isPresent),:);
end


%% Heatmaps to compare essential genes among tasks and models

% show essential genes for each model for a single specified task
figure
specTask = "Growth on Ham's media (biomass production)";
eMat_temp = squeeze(eMat(:,ismember(taskNames,specTask),:));
X = 1;  % only keep genes essential for the task in at least X models
keepGene = (sum(eMat_temp ~= 0,2) >= X);
eMat_temp = eMat_temp(keepGene,:);
genes2 = genes(keepGene);
keepGene = (sum(eMat_temp ~= 0,2)  < 27 ); % to remove common essential genes of all the models 
genHeatMap(eMat_temp(keepGene,:),modelIDs,genes2(keepGene),'both','euclidean',custom_cmap('blue'));
title(['Essential genes for task: ',specTask]);
%custom_cmap('preview');

%these results are saved as a txt file to be imprted in R
%   'EssGenesMatrix.txt'
% EssGenesMat = eMat_temp(keepGene,:);
% commonEssGeneNames = genes2(keepGene);
% modelIDs =modelIDs;
% 
%     [grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules,'Name','ENSG');
%     ihuman.grRules = grRules;
%     ihuman.genes = genes;
%     ihuman.rxnGeneMat = rxnGeneMat;
% indx =find(contains(ihuman.grRules,'TDO2'))
% for i=1:length(commonEssGeneNames)
%     indx=(find(contains(ihuman.grRules,commonEssGeneNames{i})));
%     ihuman.subSystems(indx)
% end
% 
%     for i=1:length(ihuman.subSystems)
%         ihuman.subSystems{i}=char(ihuman.subSystems{i});
%         
%     end
% 
% ihuman.subSystems(getIndexes(ihuman,commonEssGeneNames,'genes'))
% ess = readtable('EssGenesMatrix2.txt')

%this part only works if you run essentiall gene analysis on more than 1 task
% show number of essential genes for each task and each model
figure
eMat_temp = squeeze(sum(eMat,1));
X = 1;  % only keep tasks with at least X essential genes in any model
keepTask = any(eMat_temp >= X,2);
p =genHeatMap(eMat_temp(keepTask,:),modelIDs,taskNames(keepTask),'both','euclidean',custom_cmap('magma'));
title({'Number of essential genes','for each task and model'});
h = colorbar; set(get(h,'label'),'string','# Essential Genes');

% show number of tasks for which each gene is essential in each model
figure
eMat_temp = squeeze(sum(eMat,2));
X = 1;  % only keep genes essential for at least X tasks in any model
keepGene = any(eMat_temp >= X,2); 
genHeatMap(eMat_temp(keepGene,:),modelIDs,genes(keepGene),'both','euclidean',custom_cmap('blue'));
title({'Number of tasks for which each','gene is essential in each model'});
h = colorbar; set(get(h,'label'),'string','# Tasks');


% show number of models for which each gene is essential in each task
figure
eMat_temp = squeeze(sum(eMat,3));
X = 1;  % only keep tasks with at least X essential genes in any model
Y = 1;  % only keep genes essential for at least Y tasks in any model
keepTask = sum(eMat_temp ~= 0,1) >= X;
keepGene = sum(eMat_temp ~= 0,2) >= Y;
genHeatMap(eMat_temp(keepGene,keepTask),taskNames(keepTask),genes(keepGene),'both','euclidean',custom_cmap('blue'));
title({'Number of models for which each','gene is essential in each task'});
h = colorbar; set(get(h,'label'),'string','# Models');


%% Heatmaps to compare different models to each other

% compare models based on the number of tasks for which each gene is essential
figure
eMat_temp = squeeze(sum(eMat,2));
model_dist = squareform(pdist(eMat_temp','euclidean'));
genHeatMap(model_dist,modelIDs,modelIDs,'both','euclidean',custom_cmap('blue'));
title({'Difference between models in the number','of tasks for which each gene is essential'});
h = colorbar; set(get(h,'label'),'string','Euclidean distance');

% compare models based on the number of genes essential for each task
figure
eMat_temp = squeeze(sum(eMat,1));
model_dist = squareform(pdist(eMat_temp','euclidean'));
genHeatMap(model_dist,modelIDs,modelIDs,'both','euclidean',custom_cmap('blue'));
title({'Difference between models in the number','of genes essential for each task'});
h = colorbar; set(get(h,'label'),'string','Euclidean distance');


% compare models based on total Hamming distance for all tasks and genes
figure
eMat_temp = reshape(eMat,size(eMat,1)*size(eMat,2),size(eMat,3));
model_dist = squareform(pdist(eMat_temp','hamming'));
genHeatMap(model_dist,modelIDs,modelIDs,'both','euclidean',custom_cmap('blue'));
title({'Overall difference between models in the','sets of essential genes for each task'});
h = colorbar; set(get(h,'label'),'string','Hamming distance');



%% Barplots to compare models to each other

% for each task, quantify the variability among models based on the total
% hamming distance over all model pairs
figure
totalDist = arrayfun(@(i) sum(sum(triu(squareform(pdist(squeeze(eMat(:,i,:))','hamming'))))),(1:numel(taskNames))');
[~,sortInd] = sort(totalDist);
barh(totalDist(sortInd));
set(gca,'YTickLabels',taskNames(sortInd),'YTick',1:numel(taskNames),'TickLength',[0 0]);
xlabel('Total Hamming distance among models');
title({'Tasks whose essential genes varied','the most among different models'});

% for each gene, quantify the variability among models based on the total
% hamming distance over all model pairs
figure
totalDist = arrayfun(@(i) sum(sum(triu(squareform(pdist(squeeze(eMat(i,:,:))','hamming'))))),(1:numel(genes))');
keepGene = find(totalDist > 0);  % only keep genes with a nonzero distance
[~,sortInd] = sort(totalDist(keepGene));
barh(totalDist(keepGene(sortInd)));
set(gca,'YTickLabels',genes(keepGene(sortInd)),'YTick',1:numel(keepGene),'TickLength',[0 0]);
xlabel('Total Hamming distance among models');
title({'Genes with the most variability in','essentiality among different models'});

mostVarEssGenes = (genes(keepGene(sortInd)));


% filePh = fopen('mostVarEssGenesExtractedfromModelAnalysis.txt','w');
% fprintf(filePh,'%s\n',mostVarEssGenes{:});
% fclose(filePh);



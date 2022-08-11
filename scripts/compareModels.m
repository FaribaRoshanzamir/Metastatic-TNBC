clear;
setRavenSolver('gurobi')
cd('/Users/fariba/Documents/Metastatic-TNBC/models');
current= pwd;

%% comparing "tINIT Models" of Metastatic, Primary and Healthy samples
% using "compareMultipleModels" function
load('simplifiedModels_Median.mat')
comparingmodels=simplifiedModels;

clear model
for i=1:length(simplifiedModels)
    modelname=simplifiedModels{i}.modelname;
    if ~strcmpi(modelname,{'GTex_HT_Brain_CerebellarHemisphere',...
            'GTex_HT_Brain_Cerebellum',...
            'GTex_HT_Brain_NucleusAccumbens_basalGanglia',...
            'GTex_HT_Brain_Putamen_basalGanglia',...
            'GTex_HT_Brain_FrontalCortex',...
            'GTex_HT_Brain_Cortex',...
            'GTex_HT_Brain_AnteriorCingulateCortex',...
            'Brain_MediansOfAll'}) 
        disp(modelname)
        model{i}= simplifiedModels{i};
       
    end
end

%find the index of them to be excluded
find(cellfun(@isempty, model))
%excluding from final array
%clear comparingmodels
%comparingmodels=simplifiedModels([1:20 23 27:28 31 33:41]);
%length(comparingmodels)
inclIndex=find (~(cellfun(@isempty, model)))
comparingmodels=simplifiedModels(inclIndex);
length(comparingmodels)

for j=1:length(comparingmodels)
    for i=1:length(comparingmodels{j}.subSystems)
        comparingmodels{j}.subSystems{i}=char(comparingmodels{j}.subSystems{i});
        
    end
end


% seperate TP, NT, GTex, GSE dataset
for i=1:length(comparingmodels)
    modelname=comparingmodels{i}.modelname;
    if contains(modelname,'TP') 
        disp(modelname)
        modelTP{i}= comparingmodels{i};
    end
    if contains(modelname,'NT') 
        disp(modelname)
        modelNT{i}= comparingmodels{i};
    end
    if contains(modelname,'GTex') 
        disp(modelname)
        modelGTex{i}= comparingmodels{i};
    end
    if contains(modelname,'TM') 
        disp(modelname)
        modelMet{i}= comparingmodels{i};
    end 
    if ~contains(modelname,{'TM','GTex','NT','TP'}) 
        modelelse{i}= comparingmodels{i};
        disp(modelname)
    end
end

NT=find (~(cellfun(@isempty, modelNT)))
TP=find (~(cellfun(@isempty, modelTP)))
GTex=find (~(cellfun(@isempty, modelGTex)))
TM=find (~(cellfun(@isempty, modelMet)))
otherGTex = find (~(cellfun(@isempty, modelelse)))
%
cellfun(@length,{TP,NT,TM,GTex,otherGTex})


%
TNBCmetastaticProj_INITmodels([1:10])=comparingmodels(TP); %TP
TNBCmetastaticProj_INITmodels([11:14])=comparingmodels(NT);%NT
TNBCmetastaticProj_INITmodels([15:20])=comparingmodels(TM);%TM
TNBCmetastaticProj_INITmodels([21:30])=comparingmodels(GTex);%GTex
TNBCmetastaticProj_INITmodels([31:33])=comparingmodels(otherGTex);%GTex
TNBCmetastaticProj_INITmodels{31}.modelname = strcat('GTex_HT_', TNBCmetastaticProj_INITmodels{31}.modelname);
TNBCmetastaticProj_INITmodels{32}.modelname = strcat('GTex_HT_', TNBCmetastaticProj_INITmodels{32}.modelname);
TNBCmetastaticProj_INITmodels{33}.modelname = strcat('GTex_HT_', TNBCmetastaticProj_INITmodels{33}.modelname);

clear model
for i=1:length(TNBCmetastaticProj_INITmodels)
    modelname=TNBCmetastaticProj_INITmodels{i}.modelname;
    disp(modelname)
    TNBCmetastaticProj_INITmodels{i}.id=strrep(modelname, '_', ' ');
    model{i} = TNBCmetastaticProj_INITmodels{i};
end
%save('TNBCmetastaticProj_INITmodels','TNBCmetastaticProj_INITmodels');


%%
load('TNBCmetastaticProj_INITmodels.mat');

for i=1:length(TNBCmetastaticProj_INITmodels)
%Add Boundary Metabolites Because Exchange metabolites should normally not be removed from the model when using checkTasks
TNBCmetastaticProj_INITmodels{1,i} = addBoundaryMets(TNBCmetastaticProj_INITmodels{1,i});
end
% PT  models
a=zeros(1,10);
a(:)=1;
% NT models
b=zeros(1,4);
b(:)=2;
% TM models
c=zeros(1,6);
c(:)=3;
% HT models
d=zeros(1,13);
d(:)=4;
cellfun(@length,{a,b,c,d})
%
%
%
%compare multiple models
groupVector= [a,b,c,d];
length(groupVector)
%compStruct=compareMultipleModels(TNBCmetastaticProj_INITmodels,[],true,groupVector) % for test

%% Structural comparison and save the output for the further analyses
load('TNBCmetastaticProj_INITmodels.mat');
% Compare GEM structures
res = compareMultipleModels(TNBCmetastaticProj_INITmodels,false,true)
%save ('/Users/fariba/Documents/Metastatic-TNBC/Results/compareMultipleModelsOutput.mat','res')
writecell(res.modelIDs,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.modelIDs.txt','Delimiter','\t');
writecell(res.reactions.IDs,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.reactions.IDs.txt','Delimiter','\t');
writematrix(res.reactions.matrix,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.reactions.matrix.txt','Delimiter','\t');
writematrix(res.structComp,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.structComp_hammingDist.txt','Delimiter','\t');
writecell(res.subsystems.ID,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.subsystems.IDs.txt','Delimiter','\t');
writematrix(res.subsystems.matrix,'/Users/fariba/Documents/Metastatic-TNBC/Results/compStruct.subsystems.matrix.txt','Delimiter','\t');

% to have metabolite inclusion matrix you need to change
% compareMultipleModels function line 160 from "field = 'rxns';" to "field = 'metNames';"
% now you have metabolite inclusion information instead of reactions
% res2 = compareMultipleModels(TNBCmetastaticProj_INITmodels,false,true)
% writecell(res2.reactions.IDs,'/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Results/compStruct.metabolites.IDs.txt','Delimiter',',');
% writematrix(res2.reactions.matrix,'/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Results/compStruct.metabolites.matrix.txt','Delimiter',',');
% 

clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);
rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 5);
% plot and label the GEMs in tSNE space
scatter(rxn2Dmap(:,1), rxn2Dmap(:,2));
hold on
text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);


%To simplify the interpretation of results, select a subset of GEMs to analyze further.

%keep models group by metastatic tumors, primary, healthy,...
useModels = {'basalTNBC TP','BRCA basalTNBC TP','BRCA NT','GTex HT Breast','GTex HT liver','LIHC NT Liver',...
    'LIHC TP Liver', 'Liver TM'};
keep = ismember(res.modelIDs, useModels);
subMat = res.subsystems.matrix(:, keep);

%Calculate the percent difference of GEM subsystem coverage
%(number of reactions in the subsystem) from the mean coverage.

subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;

%Visualize the difference in subsystem coverage with a clustergram, including 
%only subsystems with at least a 25% difference in one or more GEMs.

% select subsystems to include in plot
inclSub = any(abs(subCoverage) > 25, 2);
subNames = res.subsystems.ID(inclSub);

% generate clustergram
cg = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, 'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels', useModels, 'ShowDendrogram', 'OFF');

%for all models
keep = ismember(res.modelIDs, res.modelIDs);
subMat = res.subsystems.matrix(:, keep);
%Calculate the percent difference of GEM subsystem coverage
%(number of reactions in the subsystem) from the mean coverage.

subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;

%Visualize the difference in subsystem coverage with a clustergram, including 
%only subsystems with at least a 25% difference in one or more GEMs.

% select subsystems to include in plot
inclSub = any(abs(subCoverage) > 25, 2);
subNames = res.subsystems.ID(inclSub);

% generate clustergram
cg = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, 'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels', res.modelIDs, 'ShowDendrogram', 'OFF');

%% Functional comparison 
% this section will take time if would be done for all 33 models

% uncomment each section to have specific comparison for each metastatic model

% useModels = {'basalTNBC TP','BRCA basalTNBC TP','BRCA NT','GTex HT Breast','GTex HT liver','LIHC NT Liver',...
%     'LIHC TP Liver', 'Liver TM'};
% keep = ismember(res.modelIDs, useModels);

useModels = {'basalTNBC TP','BRCA basalTNBC TP','BRCA NT','GTex HT Breast','Lung TM',...
    'GTex HT lung','LUSC NT Lung','LUSC TP Lung', 'LUAD NT Lung','LUAD TP Lung' };
keep = ismember(res.modelIDs, useModels);

% Task file including 57 tasks
compStruct=compareMultipleModels(TNBCmetastaticProj_INITmodels,[],true,groupVector,true,'metabolicTasks_Essential.xlsx')
compStruct2=compareMultipleModels(TNBCmetastaticProj_INITmodels,[],true,groupVector,true,'metabolicTasks_Full.xlsx')

% Task file including 57 essential tasks
% taskfile = [current '/Users/fariba/Documents/Metastatic-TNBC/data/metabolicTasks_Essential.xlsx'];

% Task file including 256 tasks
taskfile = [current '/Users/fariba/Documents/Metastatic-TNBC/data/metabolicTasks_Full.xlsx'];

res_func = compareMultipleModels(TNBCmetastaticProj_INITmodels(keep), false, false, [], true, taskfile);

%Identify which tasks differed among the GEMs (i.e., not all passed or all failed).

isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2);
diffTasks = res_func.funcComp.tasks(isDiff)

%Generate a scatter plot to visualize GEM performance on the subset of tasks
%that differed.

% visualize the matrix
spy(res_func.funcComp.matrix(isDiff,:), 30);
% apply some formatting changes
set(gca, 'XTick', 1:numel(useModels), 'XTickLabel', useModels, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');


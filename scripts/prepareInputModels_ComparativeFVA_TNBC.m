%prepareInputModels_comparativeFVA_TNBC.m

% use this script to organaize the reconstructed ecModels to have them in
% the same order as each other
%%
setRavenSolver('gurobi')
cd('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median');
current= pwd;
%%
%=================**** Load GEMs ****=========================================================================================================================
load('simplifiedModels_Median.mat')
comparingmodels=simplifiedModels;

clear model
for i=1:length(simplifiedModels)
    modelname=extractAfter(simplifiedModels{i}.description,["for "]);
    if ~strcmpi(modelname,{'GTex_Brain_CerebellarHemisphere',...
            'GTex_Brain_Cerebellum',...
            'GTex_Brain_NucleusAccumbens_basalGanglia',...
            'GTex_Brain_Putamen_basalGanglia',...
            'GTex_Brain_FrontalCortex',...
            'GTex_Brain_Cortex',...
            'GTex_Brain_AnteriorCingulateCortex',...
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

for i=1:length(comparingmodels)
    modelname=extractAfter(comparingmodels{i}.description,["for "]);
    comparingmodels{i}.id=modelname;
    disp(modelname)
    % seperate TP, NT, GTex, GSE dataset
    %modelname=comparingmodels{i}.id;
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
    if contains(modelname,'Met') 
        disp(modelname)
        modelMet{i}= comparingmodels{i};
    end 
    if ~contains(modelname,{'Met','GTex','NT','TP'}) 
        modelelse{i}= comparingmodels{i};
        disp(modelname)
    end
end

modelelse
NT=find (~(cellfun(@isempty, modelNT)))
TP=find (~(cellfun(@isempty, modelTP)))
GTex=find (~(cellfun(@isempty, modelGTex)))
TM=find (~(cellfun(@isempty, modelMet)))
otherGTex = find (~(cellfun(@isempty, modelelse)))
%
cellfun(@length,{TP,NT,TM,GTex,otherGTex})
%
TNBCmetastaticProj_INITmodels([1:length(TP)])=comparingmodels(TP); %TP
TNBCmetastaticProj_INITmodels([(length(TP)+1):(length(TP)+length(TM))])=comparingmodels(TM);%TM
TNBCmetastaticProj_INITmodels([(length(TP)+length(TM)+1):(length(TP)+length(TM)+length(NT))])=comparingmodels(NT);%NT
TNBCmetastaticProj_INITmodels([(length(TP)+length(TM)+length(NT)+1):(length(TP)+length(TM)+length(NT)+length(GTex))])=comparingmodels(GTex);%GTex
TNBCmetastaticProj_INITmodels([(length(TP)+length(TM)+length(NT)+length(GTex)+1):(length(TP)+length(TM)+length(NT)+length(GTex)+length(otherGTex))])=comparingmodels(otherGTex);%GTex
TNBCmetastaticProj_INITmodels{31}.id = strcat('GTex_', TNBCmetastaticProj_INITmodels{31}.id);
TNBCmetastaticProj_INITmodels{32}.id = strcat('GTex_', TNBCmetastaticProj_INITmodels{32}.id);
TNBCmetastaticProj_INITmodels{33}.id = strcat('GTex_', TNBCmetastaticProj_INITmodels{33}.id);

%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/tmpFiles/TNBCmetastaticProj_INITmodels_toRunFVA','TNBCmetastaticProj_INITmodels');

% clear model
% for i=1:length(TNBCmetastaticProj_INITmodels)
%     modelname=TNBCmetastaticProj_INITmodels{i}.id;
%     disp(modelname)
%     TNBCmetastaticProj_INITmodels{i}.id=strrep(modelname, '_', ' ');
%     model{i} = TNBCmetastaticProj_INITmodels{i};
% end


%=================**** Load ecGEMs ****=========================================================================================================================
clear 
load ('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/Final_ecModels/TNBC_ecModels_prT_array_1.mat')

for i=1:length(TNBC_ecModels_prT)
    modelname=TNBC_ecModels_prT{i}.name;
    if ~strcmpi(modelname,{'GTex_Brain_CerebellarHemisphere'
            ,'GTex_Brain_Cerebellum'
            ,'GTex_Brain_NucleusAccumbens_basalGanglia'
            ,'GTex_Brain_Putamen_basalGanglia'
            ,'GTex_Brain_FrontalCortex'
            ,'GTex_Brain_Cortex'
            ,'GTex_Brain_AnteriorCingulateCortex'
            ,'Brain_MediansOfAll'}) 
        disp(modelname)
        model{i}= TNBC_ecModels_prT{i};
       
    end
end

%find the index of them to be excluded
find(cellfun(@isempty, model))
%excluding from final array
%clear comparingmodels
%comparingmodels=simplifiedModels([1:20 23 27:28 31 33:41]);
%length(comparingmodels)
inclIndex=find (~(cellfun(@isempty, model)))
comparingmodels=TNBC_ecModels_prT(inclIndex);
length(comparingmodels)

for j=1:length(comparingmodels)
    for i=1:length(comparingmodels{j}.subSystems)
        comparingmodels{j}.subSystems{i}=char(comparingmodels{j}.subSystems{i});
        
    end
end

% seperate TP, NT, GTex, GSE dataset
for i=1:length(comparingmodels)
    modelname=comparingmodels{i}.id;
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
    if contains(modelname,'Met') 
        disp(modelname)
        modelMet{i}= comparingmodels{i};
    end 
    if ~contains(modelname,{'Met','GTex','NT','TP'}) 
        modelelse{i}= comparingmodels{i};
        disp(modelname)
    end
end
modelelse
NT=find (~(cellfun(@isempty, modelNT)))
TP=find (~(cellfun(@isempty, modelTP)))
GTex=find (~(cellfun(@isempty, modelGTex)))
TM=find (~(cellfun(@isempty, modelMet)))
otherGTex = find (~(cellfun(@isempty, modelelse)))
%
cellfun(@length,{TP,NT,TM,GTex,otherGTex})
%
TNBCmetastaticProj_ecModels([1:length(TP)])=comparingmodels(TP); %TP
TNBCmetastaticProj_ecModels([(length(TP)+1):(length(TP)+length(TM))])=comparingmodels(TM);%TM
TNBCmetastaticProj_ecModels([(length(TP)+length(TM)+1):(length(TP)+length(TM)+length(NT))])=comparingmodels(NT);%NT
TNBCmetastaticProj_ecModels([(length(TP)+length(TM)+length(NT)+1):(length(TP)+length(TM)+length(NT)+length(GTex))])=comparingmodels(GTex);%GTex
TNBCmetastaticProj_ecModels([(length(TP)+length(TM)+length(NT)+length(GTex)+1):(length(TP)+length(TM)+length(NT)+length(GTex)+length(otherGTex))])=comparingmodels(otherGTex);%GTex
TNBCmetastaticProj_ecModels{31}.name = strcat('GTex_', TNBCmetastaticProj_ecModels{31}.name);
TNBCmetastaticProj_ecModels{32}.name = strcat('GTex_', TNBCmetastaticProj_ecModels{32}.name);
TNBCmetastaticProj_ecModels{33}.name = strcat('GTex_', TNBCmetastaticProj_ecModels{33}.name);
%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/tmpFiles/TNBCmetastaticProj_ecModels_toRunFVA','TNBCmetastaticProj_ecModels');

% clear model
% for i=1:length(TNBCmetastaticProj_ecModels)
%     modelname=TNBCmetastaticProj_ecModels{i}.id;
%     disp(modelname)
%     TNBCmetastaticProj_ecModels{i}.id=strrep(modelname, '_', ' ');
%     model{i} = TNBCmetastaticProj_ecModels{i};
% end

%% Load GEM and ecGEM and select specific subsets
clear
%load GEMs
load('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/tmpFiles/TNBCmetastaticProj_INITmodels_toRunFVA.mat')
%load ecGEMs
load('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/tmpFiles/TNBCmetastaticProj_ecModels_toRunFVA.mat')
%=========================================================================================================================
modelIDs.tINIT = {};
for i = 1:numel(TNBCmetastaticProj_INITmodels)
    if ~ischar(TNBCmetastaticProj_INITmodels{i}.id)  % to deal with non-character IDs (cells, strings, etc)
        modelIDs.tINIT{i,1} = TNBCmetastaticProj_INITmodels{i}.id{1};
    else
        modelIDs.tINIT{i,1} = TNBCmetastaticProj_INITmodels{i}.id;
    end
end

modelIDs.gecko = {};
for i = 1:numel(TNBCmetastaticProj_ecModels)
    if ~ischar(TNBCmetastaticProj_ecModels{i}.name)  % to deal with non-character IDs (cells, strings, etc)
        modelIDs.gecko{i,1} = TNBCmetastaticProj_ecModels{i}.name{1};
    else
        modelIDs.gecko{i,1} = TNBCmetastaticProj_ecModels{i}.name;
    end
end
[tf, idx] = ismember(modelIDs.tINIT,modelIDs.gecko);
%Change order of ecModels in TNBCmetastaticProj_ecModels to be similar to
%tINITmodels in  TNBCmetastaticProj_INITmodels
TNBCmetastaticProj_ecModels_fixed =TNBCmetastaticProj_ecModels(idx);

%Select only primary and metastatic tumors
selected_ecModels =TNBCmetastaticProj_ecModels_fixed([1:16]);
selected_tINITmodels =TNBCmetastaticProj_INITmodels([1:16]);
clear comparingmodelsReArranged_fixed TNBCmetastaticProj_INITmodels comparingmodelsReArranged

cd (current)
save ('selected_ecModels.mat', 'selected_ecModels')
save ('selected_tINITmodels.mat', 'selected_tINITmodels')
%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/selected_ecModels.mat','selected_ecModels');
%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/inputFiles/selected_tINITmodels.mat','selected_tINITmodels');

%save to the compFVA folder to be used  on  Hebbe cluster
%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/compFVA_hebbe/selected_ecModels.mat','selected_ecModels');
%save('/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/run_ComparativeFVAonHebbe/compFVA_hebbe/selected_tINITmodels.mat','selected_tINITmodels');


%gen_TNBCprj_ecModels.m
%cd '/Users/fariba/Documents/GitHub'
% system('git clone https://github.com/SysBioChalmers/EnzymeConstrained_humanModels.git');
% cd EnzymeConstrained_humanModels
% system('git checkout d000641')
% cd ..
% system('git clone https://github.com/SysBioChalmers/cobratoolbox.git');
% cd cobratoolbox
% system('git checkout ...') %needed to be checked on cobratoolbox github
% repository devel branch
% cd ..
% cd /Users/fariba/Documents/GitHub
% system('git clone https://github.com/SysBioChalmers/RAVEN.git');
% cd RAVEN
% system('git checkout cc3f6e3')
% cd ..
% system('git clone https://github.com/SysBioChalmers/Human-GEM.git');
% cd Human-GEM
% system('git checkout 74de28c')
%cd ..
%**NOTICE** for running on the cluster, copy these specific versions of the used 
%packages which are in /Users/fariba/Documents/GitHub in local system to the cluster


%Add them to path
addpath(genpath('/Users/faribar/Documents/GitHub/RAVEN'))
addpath(genpath('/Users/faribar/Documents/GitHub/human-GEM'))
addpath(genpath('/Users/faribar/Documents/GitHub/cobratoolbox'))
addpath(genpath('/Users/faribar/Documents/GitHub/EnzymeConstrained_humanModels'))

%Clone GECKO and substitute human-specific scripts
% need to have GECKO_humanFiles folder in the GitHub folder
system('git clone https://github.com/SysBioChalmers/GECKO.git --branch v1.3.5')
%Replace scripts in GECKO:
cd ('/Users/fariba/Documents/GitHub')
fileNames = dir('GECKO_humanFiles');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
        fullName   = ['GECKO_humanFiles/' fileName];
        GECKOpath = dir(['GECKO/**/' fileName]);
        GECKOpath = GECKOpath.folder;
        copyfile(fullName,GECKOpath)
    end
end

%% prepare GECKO inputs
% simplify models
cd '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median'
current=pwd;
files = dir(fullfile(current, '*.mat'));
L = length(files);
mkdir simplifiedModels
for i=1:L
clear init_model;
    file=files(i).name;
    filepath = fullfile(current, file);
    load (filepath);
    modeldescription    = init_model.description; 
    modelname           = extractAfter(modeldescription,["for "]);
    AllinitModels{i}    = init_model;     
    simplifiedmodel= simplifyModel(init_model);
    [simplifiedmodel, deletedReactions, deletedMetabolites]=simplifyModel(init_model);
    % to use in GECKO pipeline model.b should be one column zero
    simplifiedmodel.b(:,2) = [];
    %%save simplified models to use by GECKO and also model comparison
    save([current '/simplifiedModels/simplifiedtINIT_model_' num2str(i) '_' modelname '.mat'],'simplifiedmodel')
    simplifiedModels{i}=simplifiedmodel;
    save('simplifiedModels_Median','simplifiedModels');
end
%save('simplifiedModelorderedBythis' ,'files')


% load tINIT models 
load ( '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median/simplifiedModels_Median.mat')
%--------------------------------------------------------------------------
%set variables to Run GECKO
org_name     = 'homo sapiens';
keggCode     = 'hsa';
GECKO_path   = '/Users/fariba/Documents/GitHub/GECKO'; % ADD GECKO PATH

L = length(simplifiedModels);
cd ('/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median') % CD TO csb_HEK293 REPOSITORY
current   = pwd;
%mkdir ecModels

%% Run GECKO
for i=1:L
clear simpModel;
    cd (current);
    simpModel = simplifiedModels{i};
    
    %------------------------- model preprocessing ------------------------
    model_modified = modelModifications(simpModel);
    modeldescription = model_modified.description;
    modelname = extractAfter(modeldescription,["for "])
    
    %------------------------- GECKO modifications ------------------------
    cd ([GECKO_path '/geckomat/get_enzyme_data'])
    model_data = getEnzymeCodes(model_modified);
    kcats = matchKcats(model_data,org_name);
    %Save ecModel data
    model_data = removeFields(model_data);
    cd ([GECKO_path '/geckomat/change_model'])
    ecModel = readKcatData(model_data,kcats);
    ecModel.name = modelname;
   
    %------------------------- Save models --------------------------------
    cd (current)
    save([current '/ecModels/' 'ecModel' '_' modelname '.mat'],'ecModel')
    %TNBC_ecModels.ecModels{i} = ecModel; 
%save ( [current '/ecModels/' 'TNBC_ecModels.mat'],'TNBC_ecModels')
end

ecDir = [current '/ecModels']
files = dir(fullfile(ecDir, '*.mat'));
L = length(files)
%mkdir arraydataFile
cd (ecDir)

for i=1:L
    file=files(i).name;
    filepath = fullfile(ecDir, file);
    load (filepath);
    modeldescription    = ecModel.description; 
    modelname           = extractAfter(modeldescription,["for "]);
    ecModel.description = [ 'Automatically generated ecModel for ' modelname ];
    ecModel.id = [modelname '_ecModel'];
    TNBC_ecModels.ecModels{i,1} = ecModel;
end
save([ecDir '/arraydataFile/' 'TNBC_ecModels_array.mat'],'TNBC_ecModels')

%% constrain using total protein 1  : this part was selected for Flux based analysis of the paper 
% use the cloned GECKO repository which GECKO_humanFiles folder from
% EnzymeConstrained_humanModels repository was copied and replaced in GECKO
% repo (lines 31-45 in this script)
load '/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median/ecModels/arraydataFile/TNBC_ecModels_array.mat'
cd '/Users/faribar/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/Final_ecModels'
current = pwd;
mkdir ecModels_ConstTotalPr1
mkdir simplified_ecModels_ConstTotalPr1

GECKO_path ='/Users/fariba/Documents/GitHub/GECKO'
for i=1:numel(TNBC_ecModels.ecModels)
    % Run the following from inside the "limit_proteins" folder in GECKO
    % (github branch: constrainEnzymesDraw )
    % Constrain the model based on RNAseq 
    % Constrain protein pool
    ecModel=TNBC_ecModels.ecModels{i,1};
    ecModel.description = [ 'Automatically generated ' ecModel.id ' constrained with total protein' ];

    cd ([GECKO_path '/geckomat'  '/limit_proteins'])
    Ptotal       = 0.593; % Human biomass total protein content [g prot/gDw]
    protCoverage = 0.5;
    sigma        = 0.5;
    [ecModel_prT,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);

    %save model
    save([current '/ecModels_ConstTotalPr1/' ecModel.id '_prT.mat'],'ecModel_prT')
    % Make simulatable
    ecModel_prT_Sim = simplifyModel(ecModel_prT,true,false,true,true);
    %save model
    save([current '/simplified_ecModels_ConstTotalPr1/'  ecModel.id '_prT_Sim.mat'],'ecModel_prT_Sim')
    % add to arrayData file
    TNBC_ecModels.ecModels_ptotConstr{i,1} = ecModel_prT; 
    TNBC_ecModels.SimpEcModels_ptotConstr{i,1} = ecModel_prT_Sim; 
    
end

TNBC_ecModels_prT = TNBC_ecModels.ecModels_ptotConstr
save ([current '/TNBC_ecModels_prT_array_1.mat'], 'TNBC_ecModels_prT');

TNBC_ecModels_prT_Sim = TNBC_ecModels.SimpEcModels_ptotConstr
save ([current '/TNBC_ecModels_prT_Sim_array_1.mat'], 'TNBC_ecModels_prT_Sim');

%% constrain using total protein 2
   % Run the following from inside the "limit_proteins" folder in GECKO
    % (github branch: constrainEnzymesDraw )

load '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/tINIT_models_Median/ecModels/arraydataFile/TNBC_ecModels_array.mat'
cd '/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/Models/Final_ecModels'
current = pwd;
mkdir ecModels_ConstTotalPr2
mkdir simplified_ecModels_ConstTotalPr2

GECKO_path ='/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/git_repos/GECKO_constrainEnzyesDraw/GECKO'
for i=1:numel(TNBC_ecModels.ecModels)
    % Run the following from inside the "limit_proteins" folder in GECKO
    % (github branch: constrainEnzymesDraw )
    % Constrain the model based on RNAseq 
    % Constrain protein pool
    ecModel=TNBC_ecModels.ecModels{i,1};
    ecModel.description = [ 'Automatically generated ' ecModel.id ' constrained with total protein' ];

    cd ([GECKO_path '/geckomat'  '/limit_proteins'])
    Ptotal       = 0.593; % Human biomass total protein content [g prot/gDw]
    total_protein_mass = Ptotal;
    non_measured = ones(length(ecModel.enzymes),1);
    ecModel_prT = constrainPool(ecModel,non_measured,total_protein_mass);
    %save model
    save([current '/ecModels_ConstTotalPr2/' ecModel.id '_prT.mat'],'ecModel_prT')
    % Make simulatable
    ecModel_prT_Sim = simplifyModel(ecModel_prT,true,false,true,true);
    %save model
    save([current '/simplified_ecModels_ConstTotalPr2/'  ecModel.id '_prT_Sim.mat'],'ecModel_prT_Sim')
    % add to arrayData file
    TNBC_ecModels.ecModels_ptotConstr{i,1} = ecModel_prT; 
    TNBC_ecModels.SimpEcModels_ptotConstr{i,1} = ecModel_prT_Sim; 
    
end
%save ([current '/TNBC_ecModels_prT_array.mat'], 'TNBC_ecModels');


TNBC_ecModels_prT = TNBC_ecModels.ecModels_ptotConstr
save ([current '/TNBC_ecModels_prT_array2.mat'], 'TNBC_ecModels_prT');

TNBC_ecModels_prT_Sim = TNBC_ecModels.SimpEcModels_ptotConstr
save ([current '/TNBC_ecModels_prT_Sim_array2.mat'], 'TNBC_ecModels_prT_Sim');


%% constrain using RNAseq data instead of proteomics data!
% Run the following from inside the "limit_proteins" folder in GECKO
% (github branch: constrainEnzymesDraw )
% Constrain the model based on RNAseq

%load geneExpression arrayData
load('/Users/fariba/Documents/BrCprojectNewdataset/MetastaticTNBC_Modeling/Data/bC2_refbGTex_nonlog_median_ArrayData.mat');
data1=arrayData{1}.levels;
gIDs=arrayData{1}.genes;

pIDs2 =ecModel_batch.enzymes(find(contains(ecModel_batch.enzGenes,gIDs)));
gIDs2 =ecModel_batch.enzGenes(find(contains(ecModel_batch.enzGenes,gIDs)));

data_temp=zeros(0,1);
for i=1:length(gIDs2)
    match=false;
    for j=1:length(arrayData{1}.genes)
        if strcmp(gIDs2{i},arrayData{1}.genes{j})&& ~match
            data_temp(i,1)=arrayData{1}.levels(j);
            %gID_temp(i,1)=arrayData{1}.genes(j);
            match=true;
        end
    end
end
% Constrain total protein pool
Ptotal       = 0.593; % Human biomass total protein content [g prot/gDw]
%Ptotal       = 0.67; % Human biomass total protein content [g prot/gDw]
protCoverage = 0.5;
sigma        = 0.5;

data_temp2 =data_temp+1
data_temp2 =log2(data_temp2)


data_temp2 =data_temp2*1000
max(data_temp2)
min(data_temp2)

[ecModel_draw,sol,MW_out,counter_out] = constrainEnzymesDraw(ecModel_D,Ptotal,sigma,[],pIDs2,data_temp2,[],[]);


% Esure the constrained model is solvable
gRate =0.1;
c_source='HMR_9034';
c_UptakeExp=1000;
sol = solveLP(ecModel_draw);



cd /Users/fariba/Documents/GitHub/GECKO_fixflexiblizeProtein_Branch/geckomat/limit_proteins
[model,enzUsages,modifications] = flexibilizeProteins(ecModel_draw,gRate,c_UptakeExp,c_source);
 solveLP(model)

cd (['/Users/faribar/Documents/GitHub/EnzymeConstrained_humanModels/ComplementaryScripts/' 'GECKO_humanFiles'])
 [ecModel_pto,enzUsages,modifications] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);
solveLP(ecModel_pto)





%% setRavenSolver('gurobi');
% use the following model to Run FVA 

load('TNBC_ecModels_prT_array1.mat', 'TNBC_ecModels_prT')
for i =1:numel(TNBC_ecModels_prT)
    clear model
    model = TNBC_ecModels_prT{i};
    solveLP(model)
    FBAsol_1_unSim.f{i,1}= ans.f
    FBAsol_1_unSim.f{i,2} =TNBC_ecModels_prT{i}.id
    model   = setHamsMedium_Jon(model,true);
    %model   = setExchangeBounds(model,mediaComps);
    solveLP(model)
    FBAsol_1_unSim.f{i,3}= ans.f
    FBAsol_1_unSim.f{i,4} =TNBC_ecModels_prT{i}.id

end

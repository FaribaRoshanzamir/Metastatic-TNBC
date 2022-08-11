% gen_TNBCprj_ecModels.m

% cd '/Users/fariba/Documents/Metastatic-TNBC/github'
% system('git clone https://github.com/SysBioChalmers/EnzymeConstrained_humanModels.git');
% cd EnzymeConstrained_humanModels
% system('git checkout d000641')
% cd ..
% system('git clone https://github.com/opencobra/cobratoolbox.git');
% cd cobratoolbox
% system('git checkout 9251835') 
% cd ..
% system('git clone https://github.com/SysBioChalmers/RAVEN.git');
% cd RAVEN
% system('git checkout cc3f6e3')
% cd ..
% system('git clone https://github.com/SysBioChalmers/Human-GEM.git');
% cd Human-GEM
% system('git checkout 74de28c')
% cd ..

%Add them to the path
addpath(genpath('/Users/fariba/Documents/Metastatic-TNBC/github/RAVEN'))
addpath(genpath('/Users/fariba/Documents/Metastatic-TNBC/github/human-GEM'))
addpath(genpath('/Users/fariba/Documents/Metastatic-TNBC/github/cobratoolbox'))
addpath(genpath('/Users/fariba/Documents/Metastatic-TNBC/github/EnzymeConstrained_humanModels'))

%Clone GECKO and substitute human-specific scripts
% need to have GECKO_humanFiles folder in the GitHub folder
system('git clone https://github.com/SysBioChalmers/GECKO.git --branch v1.3.5')

% put "GECKO_humanFiles" folder in the local GitHub folder on your
% computer. GECKO_humanFiles is available in .../GitHub/Metastatic-TNBC/data 
source = '/Users/fariba/Documents/Metastatic-TNBC/data/GECKO_HumanFiles';
mkdir /Users/fariba/Documents/Metastatic-TNBC/github/GECKO_HumanFiles
destination = '/Users/fariba/Documents/Metastatic-TNBC/github/GECKO_HumanFiles';
copyfile(source,destination)


%Replace scripts in GECKO:
cd ('/Users/fariba/Documents/Metastatic-TNBC/github')
fileNames = dir('GECKO_HumanFiles');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~strcmp(fileName,'.') && ~strcmp(fileName,'..') && ~strcmp(fileName,'.DS_Store')
        fullName   = ['GECKO_HumanFiles/' fileName];
        GECKOpath = dir(['GECKO/**/' fileName]);
        GECKOpath = GECKOpath.folder;
        copyfile(fullName,GECKOpath)
    end
end

%% prepare GECKO inputs
% simplify models
cd '/Users/fariba/Documents/Metastatic-TNBC/data/initial_Models/tINIT_models_Median'
current=pwd;
files = dir(fullfile(current, '*.mat'));
L = length(files);
for i=1:L
clear init_model;
    file=files(i).name;
    filepath = fullfile(current, file);
    load (filepath);
    modeldescription    = init_model.description; 
    modelname           = extractAfter(modeldescription,["for "]);
    %AllinitModels{i}    = init_model;     
    simplifiedmodel= simplifyModel(init_model);
    [simplifiedmodel, deletedReactions, deletedMetabolites]=simplifyModel(init_model);
    % to use in GECKO pipeline model.b should be one column zero
    simplifiedmodel.b(:,2) = [];
    % changing certain part of model names to be consistent with their ids in the
    % paper
    if contains(modelname,'Met')
    % change name ids for metastatic tumors  'Met' -> '_TM'
        modelname = regexprep(modelname,'Met','_TM');
    end 
    if contains(modelname,'GTex_')
    % change name ids for metastatic tumors  'Met' -> '_TM'
        modelname = regexprep(modelname,'GTex_','GTex_HT_');
    end 
    simplifiedmodel.modelname= modelname;
    %%save simplified models to use by GECKO and also model comparison
    save(['../../../models/model_' modelname '.mat'],'simplifiedmodel')
    simplifiedModels{i}=simplifiedmodel;
    save('../../../models/simplifiedModels_Median','simplifiedModels');
end

% changing certain part of model names to be consistent with their ids in the
% paper
% cd '/Users/fariba/Documents/Metastatic-TNBC/models';
% current=pwd;
% 
% files = dir(fullfile(current, '*.mat'));
% L = length(files);
% for i=1:L
%     oldFname=files(i).name
%     filepath_old = fullfile(current, oldFname)
%     if contains(oldFname,'Met')
%         % change name ids for metastatic tumors  'Met' -> '_TM'
%         newFname = regexprep(oldFname,'Met','_TM');
%         filepath_new = fullfile(current, newFname) 
%         movefile (filepath_old,filepath_new)
%     end
%     if contains(oldFname,'GTex_')
%     % change name ids for metastatic tumors  'Met' -> '_TM'
%         newFname = regexprep(oldFname,'GTex_','GTex_HT_');
%         filepath_new = fullfile(current, newFname)
%         movefile (filepath_old,filepath_new);
% 
%     end
% end


% load models 
load ( '/Users/fariba/Documents/Metastatic-TNBC/models/simplifiedModels_Median.mat')
%--------------------------------------------------------------------------
%set variables to Run GECKO
org_name     = 'homo sapiens';
keggCode     = 'hsa';
GECKO_path   = '/Users/fariba/Documents/Metastatic-TNBC/github/GECKO'; % ADD GECKO PATH

L = length(simplifiedModels);
cd ('/Users/fariba/Documents/Metastatic-TNBC') 
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
    modelname = simpModel.modelname
    
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
end

% Keep all the ecModels in an arraydata file
ecDir = [current '/ecModels']
files = dir(fullfile(ecDir, '*.mat'));
L = length(files)
cd (ecDir)
%mkdir arraydataFile
for i=1:L
clear ecModel
    file=files(i).name;
    filepath = fullfile(ecDir, file);
    load (filepath);
    modelname           = ecModel.name;
    if contains(modelname,'Met')
    % change name ids for metastatic tumors  'Met' -> '_TM'
        modelname = regexprep(modelname,'Met','_TM');
    end 
    if contains(modelname,'GTex_')
    % change name ids for metastatic tumors  'Met' -> '_TM'
        modelname = regexprep(modelname,'GTex_','GTex_HT_');
    end 
    ecModel.name= modelname;
    ecModel.description = [ 'Automatically generated ecModel for ' modelname ];
    ecModel.id = [modelname '_ecModel'];
    TNBC_ecModels.ecModels{i,1} = ecModel;
end
save([ecDir '/arraydataFile/' 'TNBC_ecModels_array.mat'],'TNBC_ecModels')

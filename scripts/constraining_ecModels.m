%% constrain using total protein  : this part was selected for Flux based analysis in the paper 
% use the cloned GECKO repository which GECKO_humanFiles folder from
% EnzymeConstrained_humanModels repository was copied and replaced in GECKO
% repo (lines 40_51 in gen_TNBCprj_ecModels)

load '/Users/fariba/Documents/Metastatic-TNBC/ecModels/arraydataFile/TNBC_ecModels_array.mat'
cd '/Users/fariba/Documents/Metastatic-TNBC/ecModels'
current = pwd;

GECKO_path ='/Users/fariba/Documents/Metastatic-TNBC/github/GECKO'

for i=1:numel(TNBC_ecModels.ecModels)
    % Run the following from inside the "limit_proteins" folder in GECKO
    % (github branch: constrainEnzymesDraw )
    % Constrain protein pool
    ecModel=TNBC_ecModels.ecModels{i,1};
    ecModel.description = [ 'Automatically generated ' ecModel.id ' constrained with total protein' ];

    cd ([GECKO_path '/geckomat'  '/limit_proteins'])
    Ptotal       = 0.593; % Human biomass total protein content [g prot/gDw]
    protCoverage = 0.5;
    sigma        = 0.5;
    [ecModel_prT,~,~] = constrainEnzymes(ecModel,Ptotal,sigma,protCoverage);

    %save modeles one by one
    %save([current '/ecModels_Constrained_TotalProtein/' ecModel.id '_prT.mat'],'ecModel_prT')

    % add to arrayData file
    TNBC_ecModels.ecModels_ptotConstr{i,1} = ecModel_prT; 
    
end

TNBC_ecModels_prT = TNBC_ecModels.ecModels_ptotConstr
save ([current '/arraydataFiles/TNBC_ecModels_constrPrT_array.mat'], 'TNBC_ecModels_prT');


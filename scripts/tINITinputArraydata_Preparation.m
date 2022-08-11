%
% FILE NAME: tINITinputArraydata_Preparation.m
%
% DATE CREATED: 2019-06-12

% load median gene expression  "QN_TPM_bC_nonlog_median" from /data  

%%%%%personalNOTE%%%%%
% it's old name was "QN_TPM_exonLength_bC2_refbGTex_nonlog_median_nonNegative_repeatedlinesbasedonENSids"
% and stored in
% /Users/fariba/Documents/BrCprojectNewdataset/AB/lastVersionGtex/tINITinputs/FinalFiles
%%%%%%%%%%%%%%%%%%%%%%

cd ('/Users/fariba/Documents/BrCprojectNewdataset/AB/lastVersionGtex/tINITinputs/FinalFiles');
clear;

%% extracting needed information from table
opts = detectImportOptions('QN_TPM_exonLength_bC2_refbGTex_nonlog_median_nonNegative_repeatedlinesbasedonENSids.csv');
preview('QN_TPM_exonLength_bC2_refbGTex_nonlog_median_nonNegative_repeatedlinesbasedonENSids.csv',opts)

WholeDataTable_Median = readtable('QN_TPM_exonLength_bC2_refbGTex_nonlog_median_nonNegative_repeatedlinesbasedonENSids.csv',opts);

%%%%%%%%%%%% Median
current = pwd;
WholeDataTable_Median2 =WholeDataTable_Median(:,2:end-1);
colnames = WholeDataTable_Median2.Properties.VariableNames;
L = length(colnames);
colnames(1)
for i=1:L
   %disp (i)
   arrayData{i}.tissues    = {char(colnames(i))};
   arrayData{i}.genes     = table2cell(WholeDataTable_Median(:,end));
   arrayData{i}.levels    = table2array(WholeDataTable_Median2(:,char(colnames(i))));
   arrayData{i}.threshold = 1;
end

%%%%%save
%save('bC2_refbGTex_nonlog_median_ArrayData','arrayData')
save (fullfile(current, 'arrayData','bC_nonlog_median_ArrayData' ),'arrayData')

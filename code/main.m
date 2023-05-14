clc;
clear;
close all;
rootPath = '/home/data/tianxiao/Field_Inhomogeneity_project/code/code_to_git';
dataPath = fullfile(rootPath,'data');
%% load data
load([dataPath,'/preparedData.mat']);
load([dataPath,'/Subspace.mat'])
load([dataPath,'/Subspace_VT.mat'])
load(fullfile(dataPath,'sliceProfile_FOV120.mat'));
SPcorrvec = fliplr(SPcorrvec);
%% set param
rmvMarginSlice = 3;
rank_B1 = 4;
rank_RP = 5;
TR = 0.0067;
TE = 0.0028;
FA = [2/180*pi,12/180*pi,14/180*pi];
%% VFA fitting
dataSize = size(brainMask);
[E1,E2] = multiFAflashT1Fitting(flash_multiFA,FA,brainMask);
T1_uncorrected = -E1.*TR*2;
PD_eff_uncorrected = E2;
fprintf('Finish fititng!!!');
%% Tissue field
GM_mask = zeros(dataSize);GM_mask(c1rmprage>255*0.95) = 1;
WM_mask = zeros(dataSize);WM_mask(c2rmprage>255*0.95) = 1;

T1_GM = 1.6;T1_WM = 0.9;
B1_GM = (T1_uncorrected./(T1_GM.*GM_mask)).^0.5;
B1_WM = (T1_uncorrected./(T1_WM.*WM_mask)).^0.5;
B1_GWM = B1_GM + B1_WM;B1_GWM(isinf(B1_GWM))=0;B1_GWM(isnan(B1_GWM))=0;
GWM_mask = zeros(dataSize);
GWM_mask(B1_GWM>0) = 1;GWM_mask = logical(GWM_mask);

%% B1 estimation
lambda_regB1_eigen = 0.001;
SPcorrB1 = repmat(reshape(SPcorrvec,1,1,size(SPcorrvec,2)),dataSize(1),dataSize(2),1);
SPcorrB1(isnan(SPcorrB1))=0;SPcorrB1(isinf(SPcorrB1))=0;

B1_subspace = zeros(dataSize(1:3));
[B1_subspaceVec,c_B1] = B1SubspaceRecon_EigenPenalty_Mask_SPcorr(B1_GWM(brainMask)',GWM_mask(brainMask)',S_B1,VT_B1,rank_B1,B1_projection,lambda_regB1_eigen,SPcorrB1(brainMask));
B1_subspace(brainMask) = B1_subspaceVec;
T1_corrected = T1_uncorrected./(B1_subspace.^2)./(SPcorrB1.^2);
%% RP estimation
PD_GM = 0.84;PD_WM = 0.68;
TE_isoin = 0.0075;TR_isoin = 0.025;
T2star_GM = 0.066;T2star_WM = 0.05;
RP_GM = (E2./(PD_GM.*exp(-TE_isoin/T2star_GM).*B1_subspace)).*GM_mask;
RP_WM = (E2./(PD_WM.*exp(-TE_isoin/T2star_WM).*B1_subspace)).*WM_mask;
RP_GWM = RP_GM + RP_WM;RP_GWM(isinf(RP_GWM))=0;RP_GWM(isnan(RP_GWM))=0;
GWM_mask = zeros(dataSize);
GWM_mask(RP_GWM>0) = 1;GWM_mask = logical(GWM_mask);

lambda_regRP_eigen = 0.001;
RP_subspace = zeros(dataSize(1:3));
[RP_subspaceVec,c_RP] = B1SubspaceRecon_EigenPenalty_Mask(RP_GWM(brainMask)',GWM_mask(brainMask)',S_RP,VT_RP,rank_RP,RP_projection,lambda_regRP_eigen);
RP_subspace(brainMask) = RP_subspaceVec;
PD_eff_corrected = PD_eff_uncorrected./(B1_subspace.*RP_subspace.*SPcorrB1);

B1_subspace(isnan(B1_subspace))=0;B1_subspace(isinf(B1_subspace))=0;B1_subspace(B1_subspace<0)=0;B1_subspace(B1_subspace>10)=0;
RP_subspace(isnan(RP_subspace))=0;RP_subspace(isinf(RP_subspace))=0;RP_subspace(RP_subspace<0)=0;RP_subspace(RP_subspace>10000)=0;
T1_uncorrected(isnan(T1_uncorrected))=0;T1_uncorrected(isinf(T1_uncorrected))=0;T1_uncorrected(T1_uncorrected<0)=0;T1_uncorrected(T1_uncorrected>10)=0;
T1_corrected(isnan(T1_corrected))=0;T1_corrected(isinf(T1_corrected))=0;T1_corrected(T1_corrected<0)=0;T1_corrected(T1_corrected>10)=0;
PD_eff_uncorrected(isnan(PD_eff_uncorrected))=0;PD_eff_uncorrected(isinf(PD_eff_uncorrected))=0;PD_eff_uncorrected(PD_eff_uncorrected<0)=0;
PD_eff_corrected(isnan(PD_eff_corrected))=0;PD_eff_corrected(isinf(PD_eff_corrected))=0;PD_eff_corrected(PD_eff_corrected<0)=0;PD_eff_corrected(PD_eff_corrected>10)=0;

%% save results
savePath = fullfile(dataPath,'Subspace_result');
if ~exist(savePath,'dir')
    mkdir(savePath);
end
save([savePath,'/Subspace_results_rank',num2str(rank_B1),'.mat'],'B1_subspace','RP_subspace', ...
    'T1_corrected','PD_eff_corrected');

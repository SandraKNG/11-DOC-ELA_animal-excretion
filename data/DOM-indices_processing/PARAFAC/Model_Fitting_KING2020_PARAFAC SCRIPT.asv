%-----Fitting new EEMs to Williams et al. 2013 7 Component Model-----
%By Sarah King, November 2020

%Requirements for EEMs:
%EEMs must be collected and processed in the same way as the original EEMs
%or this model is not appropriate. 

%Collected on the same Varian Cary Eclipse fluorometer
%Same excitation and emission wavelength ranges and bandwidth
%   Emission 300-600 nm, 2 nm intervals
%   Excitation 250-600 nm, 5 nm intervals
%   Note: Our current scan method measures a larger range of em and ex
%   wavelengths but at the same intervals; this is why we trim the EEMs
%   (step 4 in this script)
%Processed according to the methods in Clay Williams' Nutrient Protocols
%   Instrument spectral correction
%   Inner Filter Effects (IFE) correction
%   MilliQ blank subtraction
%   Raman normalization    

%You will need the model file (11Jan2010Model7Ex250Em300AllData.mat), the
%DOMFluor toolbox and the FittingbyYY file.

%OTHER HELPFUL INFO
%To produce plots of the seven components (do step 1 first to load model):
ComponentEEM(AnalysisData,7)

%For an explanation of the modelling process and how to evaluate the
%residuals for good model fit, read
%Stedmon and Bro, 2008. Characterizing dissolved organic matter
%   fluorescence with parallel factor analysis: a tutorial. Limnology and
%   Oceanography: Methods 6: 572-579.
%   (A detailed tutorial is also available as appendix 1)

%Step 1; Load the existing PARAFAC model. Delete B, C, FMax, X and XBackup
%        Edit filepath
load 'C:\Users\sandraklemet\Documents\PARAFAC\PARAFAC\11Jan2010Model7Ex250Em300AllData.mat'
clear ('B', 'C', 'FMax')
clear ('X', 'XBackup')
clear ('Test2')

%Step 2a; Change to directory with input data (corrected EEMs as .csv, labelled PARAFAC_0001, PARAFAC_0002,..., PARAFAC_000n), load new EEMs as X 
X = combineSampleFiles('*.csv');

%Step 2b; Change directory to DOMFluor and define parameters for structure
%         Edit Em, Ex, nEm, nEx, nSample according to your new data set 
Em=(270:2:600)';%Emission range and bandwidth (min wavelength, bandwidth, max wavelength)
Ex=(230:5:500)';%Excitation range and bandwidth (min wavelength, bandwidth, max wavelength)
nEm = 166 %Number of Em
nEx = 55 %Number of Ex
nSample = 98 %Number of samples
XBackup = X; %Just Backup for original data
FittingData = struct('Ex',Ex,'Em',Em,'X',X,'nEx',nEx,'nEm',nEm,'nSample',nSample)
clear ('Em', 'Ex', 'X', 'XBackup', 'nEm', 'nEx', 'nSample')

%Step 3; Cutting the region of the spectra influenced by scatter peaks
[CutFittingData]=EEMCut(FittingData,20,20,NaN,NaN,'No')

%Step 4; Cut EEMs so they have Em = 300:2:600 and Ex = 250:5:500
[CutFittingData]=RemoveOutliers(CutFittingData,[],[1:15],[1:4]) %(Data,[outsamples],[?ex?],[?em?])
PlotEEMby4(1,CutFittingData,'RU')

%Step 5; Change directory to nway and fit components to each EEM 
%        Edit number of samples
[A,B,C]=fac2let(AnalysisData.Model7);
A=rand(98,7);%(number of samples, components)
OldLoad={A;B;C};
OldLoad=OldLoad';
clear A B C;

%Step 6; Change directory to nway and run (make sure FittingbyYY is in there)
[FitModel,Iter,Err]=FittingbyYY(CutFittingData.X,7,[1e-6],[2 2 2],OldLoad,[0 1 1]);
[Fitted]=FittingbyYY(CutFittingData.X,7,[1e-6],[2 2 2],OldLoad,[0 1 1]);
eval(['CutFittingData.Model',int2str(7),'=FitModel',]);

%Step 7; Check the residuals to ensure good fit.
EvalModel(CutFittingData,7)

%Step 8; Export the PARAFAC result
%        Edit filepath
[FMax,B,C]=ModelOut(CutFittingData,7,'C:\Users\sandraklemet\Documents\PARAFAC\PARAFAC\Results.xlsx');

cd 'C:\Documents and Settings\xxxx\MATLAB';


clear; 
clc;  

%Read data files in
OriginalData.Ex = csvread('ex.csv',1); 
OriginalData.Em = csvread('em.csv',1); 

%fl data, (columns are excitation data, rows emission x sample)
OriginalData.X= [csvread('01f_ru.csv',1)]; 

%key data (for sample id, if you want)
OriginalData.key=[xlsread('01key.csv',1)];
              
              
OriginalData.nEx=(size(OriginalData.Ex,1)); %identifys the number of Excitation wavelengths
OriginalData.nEm=(size(OriginalData.Em,1)); %identifys the number of Emission wavelengths
OriginalData.nSample=(size(OriginalData.X,1)); OriginalData.nSample=OriginalData.nSample/OriginalData.nEm; %identifys the number of samples


%reorganises the data
OriginalData.X=permute((reshape(OriginalData.X',OriginalData.nEx,OriginalData.nEm,OriginalData.nSample)),[3 2 1]);





%plots EEMs of the data with a 0.2s pause between plots.
PlotEEMBy4(1:OriginalData.nSample,OriginalData,'');
PlotSurfBy4(100:OriginalData.nSample,OriginalData,'');

%deletes unwanted items from workspace
clear i;
%saves workspace
save mydata
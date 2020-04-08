%tightCollimated
%This file is used to load the data from running "multirun_dens_gc.sh"

%Initialization
clear; close all; clc;
addpath ~/Desktop/codes/beamLaser_Proj/tightCollimated/;
addpath ~/Desktop/codes/beamLaser_Proj/common/;

%get parameters from the input.txt in the controlType directory
prompt = "What is the controlType?\n"; 
controlType = input(prompt,'s');
cd(controlType);
getParam_multiRun_tightCollimated;%This will give us the fast value

%bash variables; 

%gc
%If only one gcList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nMaxGc = 20;%
initGc = 0.001;%
intervalGc = 0.002;
gcList = initGc+(0:nMaxGc-1)*intervalGc;%dim = 1*nMaxGc

% %If more than one gcList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %gcList1 = 0.5 ~ 1.4;
% nMaxGc1 = 40;%
% initTau1 = 0.1;%
% intervalGc1 = 0.1;
% gcList1 = (initTau1+(0:nMaxGc1-1)*intervalGc1);
% %gcList = 1:nMaxGc;%in the unit of tau, if initGc = intervalGc;
% %gcList2 = 1.0 ~ 10.5
% nMaxGc2 = 0;%20/0
% initTau2 = 1.0;
% intervalGc2 = 0.5;
% gcList2 = (initTau2+(0:nMaxGc2-1)*intervalGc2);
% %for convenience
% nMaxGc = nMaxGc1+nMaxGc2;
% gcList = [gcList1,gcList2];

%dens
%If only one densList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nMaxDens = 1;%
initDens = 1000;
intervalDens = 100;
densList =initDens+(0:nMaxDens-1)*intervalDens;%dim = 1*nMaxDens

%If more than one densList%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nMaxDens1 = 40;%
% initDens1 = 100;
% intervalDens1 = 100;
% densList1 =initDens1+(0:nMaxDens1-1)*intervalDens1;%dim = 1*nMaxDens1
% nMaxDens2 = 10;
% initDens2 = 4500;
% intervalDens2 = 500;
% densList2 =initDens2+(0:nMaxDens2-1)*intervalDens2;%dim = 1*nMaxDens2
% % For convenience
% nMaxDens = nMaxDens1+nMaxDens2;
% densList = [densList1, densList2];
% %define filenames
% filenameList1 = repmat(" ", 1, nMaxDens1);
% for j = 1:nMaxDens1
%     dens1 = initDens1+intervalDens1*(j-1);
%     filenameList1(1,j) = ['pois0_dt0.01_dZ0_dPz0_tau0.50_nBin20_dens',...
%         num2str(dens1),'_gc0.1'];
% end
% filenameList2 = repmat(" ", 1, nMaxDens2);
% for j = 1:nMaxDens2
%     dens2 = initDens2+intervalDens2*(j-1);
%     filenameList2(1,j) = ['pois0_dt0.01_dZ0_dPz0_tau0.50_nBin20_dens',...
%         num2str(dens2),'_gc0.1'];
% end
% filenameList = [filenameList1, filenameList2];

%define the filenames
filenameList = repmat(" ", nMaxGc, nMaxDens);
for j = 1:nMaxDens
    dens = densList(j);
    for i = 1:nMaxGc
        gc = gcList(i);
        filenameList(i,j) = ['dens', num2str(dens),'_gc',num2str(gc,'%.3f')];
    end
end

%define empty data structure for variables
nAtom = zeros(nMaxGc, nMaxDens, nStore);
intensity = zeros(nMaxGc, nMaxDens, nStore);
% inversionAve = zeros(nMaxGc,nMaxDens,nStore);
% if fast == 0
%     spinSpinCorAve = zeros(nMaxGc,nMaxDens,nStore);
% end
% szFinal = zeros(nMaxGc,nMaxDens,nTimeStep);
JxMatrix = zeros(nMaxGc, nMaxDens, nTrajectory, nStore);
JyMatrix = zeros(nMaxGc, nMaxDens, nTrajectory, nStore);
% szMatrix = zeros(nMaxGc,nMaxDens, nBin,nStore);

%load all the data needed
timeRequired = 3*nMaxGc*nMaxDens*(nStore/1000)*(nTrajectory/1000);
printWords1 = ['Begin loading. Time required: ' , num2str(timeRequired), ' seconds.'];
disp(printWords1);
tic;

for j = 1:nMaxDens
    for i = 1:nMaxGc
        %define the name of the directory
        cd(filenameList(i,j));
        %Load data.
        nAtom(i,j,:) = load('nAtom.dat');
        intensity(i,j,:) = load('intensity.dat');
%         inversionAve(i,j,:) = load('inversionAve.dat');
%         if fast == 0
%             spinSpinCorAve(i,j,:) = load('spinSpinCorAve_re.dat');
%         end
        %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
        %sxFinal = load('sxFinal.dat');
        %syFinal = load('syFinal.dat');
%         szFinal(i,j,:) = load('szFinal.dat');
        JxMatrix(i,j,:,:) = load('JxMatrix.dat');
        JyMatrix(i,j,:,:) = load('JyMatrix.dat');
%         szMatrix(i,j,:,:) = load('szMatrix.dat');
        %spinSpinCor_re = load('spinSpinCor_re.dat');
        %spinSpinCor_im = load('spinSpinCor_im.dat');  
        cd ..;
            %Print out info
        percentage = (i+(j-1)*nMaxGc)/nMaxGc/nMaxDens*100;
        printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
        disp(printWords2);
    end
end
% cd ..;
toc;
%     for i = 1:nMaxGc2
%         %define the name of the directory
%         tau2 = initTau2+intervalGc2*(i-1);
%         nAtomAve = initDens+intervalDens*(j-1);
%         filename = ['pois_tau', num2str(tau2,'%.1f'), '_g2_k40_N', num2str(nAtomAve)];
%         cd(filename);
%         %Load data.
%         nAtom(i+nMaxGc1,j,:) = load('nAtom.dat');
%         intensity(i+nMaxGc1,j,:) = load('intensity.dat');
%         inversionAve(i+nMaxGc1,j,:) = load('inversionAve.dat');
%         spinSpinCorAve(i+nMaxGc1,j,:) = load('spinSpinCorAve_re.dat');
%         %spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
%         %sxFinal = load('sxFinal.dat');
%         %syFinal = load('syFinal.dat');
%         szFinal(i+nMaxGc1,j,:) = load('szFinal.dat');
%         qMatrix(i+nMaxGc1,j,:,:) = load('qMatrix.dat');
%         pMatrix(i+nMaxGc1,j,:,:) = load('pMatrix.dat');
%         %spinSpinCor_re = load('spinSpinCor_re.dat');
%         %spinSpinCor_im = load('spinSpinCor_im.dat');  
%         cd ..;
%             %Print out info
%         percentage = (i+nMaxGc1+(j-1)*nMaxGc)/nMaxGc/nMaxDens*100;
%         printWords2 = ['Finish loading ', num2str(percentage,'%.1f'),'%...'];
%         disp(printWords2);
%     end
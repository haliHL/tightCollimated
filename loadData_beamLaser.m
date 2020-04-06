%cNumber_old

%Initialization
clear; close all; 
%clc;

%Get parameters
getParam_beamLaser;
projectdir1 = controlType;
cd(projectdir1);
projectdir2 = name;
cd(projectdir2);

%Load data.
fprintf('Loading data... This may take several minutes.\n')
nAtom = load('nAtom.dat');
intensity = load('intensity.dat');
% inversionAve = load('inversionAve.dat');
JxMatrix = load('JxMatrix.dat');
JyMatrix = load('JyMatrix.dat');

% jzMatrix_true = 0;
% if isfile('JzMatrix.dat') && isempty(load('JzMatrix.dat'))==false
%     jzMatrix_true = 1;
%     JzMatrix = load('JzMatrix.dat');
% end

% sFinalMatrix_true = 0;
% if isfile('sxFinalMatrix.dat') && isempty(load('sxFinalMatrix.dat'))==false
%     sFinalMatrix_true = 1;
%     sxFinalMatrix = load('sxFinalMatrix.dat');    
%     syFinalMatrix = load('syFinalMatrix.dat');    
%     szFinalMatrix = load('szFinalMatrix.dat');    
% end

if fast == 0
%     spinSpinCorAve_re = load('spinSpinCorAve_re.dat');
    % spinSpinCorAve_im = load('spinSpinCorAve_im.dat');
%     spinSpinCor_re = load('spinSpinCor_re.dat');
%     spinSpinCor_im = load('spinSpinCor_im.dat');
    sxMatrix = load('sxMatrix.dat');
    syMatrix = load('syMatrix.dat');
    szMatrix = load('szMatrix.dat');
end

%Add the previous directory into path
addpath ~/Desktop/codes/beamLaser_Proj/cNumber_old/;
addpath ~/Desktop/codes/beamLaser_Proj/common/;
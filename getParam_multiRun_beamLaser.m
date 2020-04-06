%cNumber_old
%Import simulation parameters from the input file in controlType dir
dataTable = readtable('inputControl.txt');
dataVariables = dataTable.Variables;
%Define all the simulation parameters one by one;
dt = str2double(dataVariables{1,2});
tmax = str2double(dataVariables{2,2});
nStore = str2double(dataVariables{3,2});
nTrajectory = str2double(dataVariables{4,2});
nBin = str2double(dataVariables{5,2});
yWall = str2double(dataVariables{6,2});
lambda = str2double(dataVariables{7,2});
deltaZ = str2double(dataVariables{8,2});
deltaPz = str2double(dataVariables{9,2});
transitTime = str2double(dataVariables{10,2});
density = str2double(dataVariables{11,2});
mAtom = str2double(dataVariables{12,2});
rabi = str2double(dataVariables{13,2});
kappa = str2double(dataVariables{14,2});
detuning = str2double(dataVariables{15,2});
T2 = str2double(dataVariables{16,2});
controlType = dataVariables{17,2};
name = dataVariables{18,2};
%pois = dataVariables{19,2};
fast = str2double(dataVariables{20,2});




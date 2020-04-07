%tightCollimated
%Import simulation parameters from the input file in controlType dir
dataTable = readtable('inputControl.txt');
dataVariables = dataTable.Variables;
%Define all the simulation parameters one by one;
tmax = str2double(dataVariables{1,2});
nStore = str2double(dataVariables{2,2});
nTrajectory = str2double(dataVariables{3,2});
nBin = str2double(dataVariables{4,2});
density = str2double(dataVariables{5,2});
gc = str2double(dataVariables{6,2});
controlType = dataVariables{7,2};
name = dataVariables{8,2};
fast = str2double(dataVariables{9,2});




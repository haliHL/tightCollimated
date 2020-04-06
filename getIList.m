%cNumber_full

%This file is to solve the transcendental equation
%   kappa*a^2=phi*sin^2(g*a*tau/2) 
%as a fuction of tauList.

%Input: tauList, dim = (1, N)
%Output: idaList, dim = (1, N)

function iList = getIList(tauList, gc, density)

iList = tauList.*0;
for i = 1:size(tauList,2)
    iList(i) = getIValue(tauList(i),gc,density);
end
    

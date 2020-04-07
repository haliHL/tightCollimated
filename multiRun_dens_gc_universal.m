%tightCollimated

%This code needs to be run after running "loadData_multiRun_dens_gc.m".
%This code treats data with a 2D UNIVERSAL parameter manner.
%% PREPARATIONS

transitTime = 1;
tStep = tmax / nStore;

%set the steadyMultiplier
steadyMultiplier = 5;%This value can be varied if needed. 5 is empirical.
t0 = steadyMultiplier * transitTime;
n0_nStore = ceil(t0 / tmax * nStore);
n0_nTimeStep = ceil(t0 / tmax * nTimeStep);
m = nStore - n0_nStore + 1;

%get ntgcMatrix
ntgcMatrix = zeros(nMaxGc, nMaxDens);
for j = 1:nMaxDens
    for i = 1:nMaxGc
        ntgcMatrix(i,j) = densList(j)*transitTime^2*gcList(i);
    end
end
%get ntgcList and ntgcOrder
ntgcSize = nMaxDens*nMaxGc;
[ntgcList,ntgcOrder] = sort(reshape(ntgcMatrix, [1, ntgcSize]));

%Define batch parameters
batchTime = min(100*transitTime, tmax-t0);%approximate coherence time
batchSize = batchTime/tStep;%want at lest 20
tBatch = linspace(0,batchTime,batchSize);
nBatch = floor(m/batchSize);%want at least 20
%% Analytical Solutions

%get analytical solution
iValueList = zeros(1, ntgcSize);
for i = 1:ntgcSize
    myfun = @(x, y) 0.5*(1-cos(sqrt(x*y)))-y;
    x = ntgcList(i);
    fun = @(y) myfun(x,y);
    y = fzero(fun, 0.5);
    
    iValueList(i) = y;
end
% get iValueMatrix
iValueMatrix = zeros(nMaxGc,nMaxDens);
for j = 1:nMaxDens
    for i = 1:nMaxGc
        iValueMatrix(i,j) = getIValue(transitTime,gcList(i),densList(j));
    end
end
% getIvalueList
iValueList = reshape(iValueMatrix./densList, [1, ntgcSize]);
iValueList = iValueList(ntgcOrder);
plot(ntgcList, iValueList);
%% intensity

%load data
%intensitySS
intensitySS = intensity(:,:,n0_nStore:nStore);
intensitySS_mean = mean(intensitySS, 3);
intensitySS_std = std(intensitySS, 0, 3);
%get intensitySSList
intensitySSList = reshape(intensitySS_mean./densList, [1, ntgcSize]);
intensitySSList = intensitySSList(ntgcOrder);

%plot data
figure(1);
hold on;
scatter(ntgcList,intensitySSList, '+');
% plot(ntgcList,iValueList,'LineWidth',1.5);
hold off;
ax = gca;
ax.FontSize = 20;
% ax.XLim = [0 100];
ax.YLim = [0 1];
xlabel('\Phi\tau^{2}\Gamma_c','FontSize', 20);
ylabel('I/\Phi','FontSize', 20);

fprintf('Press enter to see intensity with different tau values specified.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot data with different gc's
%{
% figure(11);
% hold on;
% scatter(ntgcMatrix(2,:),intensitySS_mean(2,:)./densList,'+'); 
% scatter(ntgcMatrix(4,:),intensitySS_mean(4,:)./densList,'o'); 
% scatter(ntgcMatrix(6,:),intensitySS_mean(6,:)./densList,'x'); 
% scatter(ntgcMatrix(8,:),intensitySS_mean(8,:)./densList,'p'); 
% % for i = 1:nMaxGc
% %     scatter(ntgcMatrix(i,:),intensitySS_mean(i,:)./densList,'+');    
% % end
% hold off;
% ax = gca;
% ax.FontSize = 20;
% ax.XLim = [0 100];
% ax.YLim = [0 1];
% legend({'\Gamma_C = 0.02'
%         '\Gamma_C = 0.04'
%         '\Gamma_C = 0.06'
%         '\Gamma_C = 0.08'});
% % legend({'\alpha = 0.2'
% %         '\alpha = 0.3'
% %         '\alpha = 0.4'
% %         '\alpha = 0.5'
% %         '\alpha = 0.6'
% %         '\alpha = 0.7'
% %         '\alpha = 0.8'
% %         '\alpha = 0.9'
% %         '\alpha = 1.0'
% %         '\alpha = 1.1'});
% xlabel('\Phi\tau^{2}\Gamma_c','FontSize', 20);
% ylabel('I/\Phi','FontSize', 20);
% 
% fprintf('Intensity shown. Press enter to see linewidth next.\n');
% pause;
%}
%% linewidth

%load data
%fMatrix
fMatrix1 = zeros(nMaxGc,nMaxDens);
fMatrix2 = zeros(nMaxGc,nMaxDens);
fMatrix3 = zeros(nMaxGc,nMaxDens);
%loop over all gc's and dens's;
%get corresponding parameters from loadData_multiRun_dens_gc.m;
for j = 1:nMaxDens
    for i = 1:nMaxGc
        %get q and p for a single simulation
        q = zeros(nTrajectory, m);
        p = zeros(nTrajectory, m);
        q(:,:) = -JyMatrix(i,j,:,n0_nStore:nStore);
        p(:,:) = JxMatrix(i,j,:,n0_nStore:nStore);
        %get non-negative part of g(1) function
        realG1Pos = zeros(1,batchSize);
        imagG1Pos = zeros(1,batchSize);
        for k=1:nBatch
            q_k = q(:,(1+batchSize*(k-1)):(k*batchSize));
            p_k = p(:,(1+batchSize*(k-1)):(k*batchSize));
            realG1Pos = realG1Pos + (q_k(:,1)'*q_k+p_k(:,1)'*p_k);
            imagG1Pos = imagG1Pos + (p_k(:,1)'*q_k-q_k(:,1)'*p_k);
        end
        %normalize g(1) with realG1Pos. Notice that the relative value matters.
        maxRealG1Pos = max(realG1Pos);
        realG1Pos = realG1Pos/maxRealG1Pos;
        imagG1Pos = imagG1Pos/maxRealG1Pos;
        
        %Method 1: exp fit
        l1 = coeffvalues(fit(tBatch',realG1Pos','exp1','startpoint',[1,-0.5])); 
        fMatrix1(i,j) = -l1(2)/pi;
        
        %Method 2: fwhm
        realG1Plot = [fliplr(realG1Pos),realG1Pos(2:end)];
        imagG1Plot = [-fliplr(imagG1Pos),imagG1Pos(2:end)];
        G1Plot = realG1Plot+1i*imagG1Plot;
        S = fft(ifftshift(G1Plot));
        timeInterval = tStep;
        df = 1/timeInterval;
        %smoothdata() equivalent to decreasing total time
        s = real(S);
        %normalize s such that the peak is at 1;
        s = fftshift(s)/max(s);
        L = size(s,2);
        f = df*(-L/2+1:L/2)/L;
        spectra = [f; s]';
        fMatrix2(i,j) = fwhm(spectra(:,1),spectra(:,2));
        fMatrix3(i,j) = fwhm2(spectra(:,1),spectra(:,2));
    end
end
%fList
fList1 = reshape(fMatrix1./gcList', [1, ntgcSize]);
fList1 = fList1(ntgcOrder);
fList2 = reshape(fMatrix2./gcList', [1, ntgcSize]);
fList2 = fList2(ntgcOrder);
fList3 = reshape(fMatrix3./gcList', [1, ntgcSize]);
fList3 = fList3(ntgcOrder);
%debug
% fList3 = reshape(linewidth_exp(2:end, 2:end)./gcList', [1, ntgcSize]);
% fList3 = fList3(ntgcOrder);
%debug

%plot data
figure(2);
hold on;
scatter(ntgcList,fList1,'+');
% scatter(ntgcList,fList2,'o');
scatter(ntgcList,fList3,'*');
%debug
% scatter(ntgcList,fList3,'x');
%debug
xlabel('\Phi\tau^{2}\Gamma_c','FontSize', 20);
ylabel('\Delta \nu/\Gamma_c','FontSize', 20);
ax = gca;
ax.FontSize = 20;
ax.YScale = 'log';
ax.XLim = [0 20];
ax.YLim = [0.01 inf];

fprintf('Press enter to see linewidth with different tau values specified.\n');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot data with different tau's
%{ 
% figure(21);
% hold on;
% scatter(ntgcMatrix(2,:),fMatrix(2,:)/gcList(2),'+');  
% scatter(ntgcMatrix(4,:),fMatrix(4,:)/gcList(4),'o'); 
% scatter(ntgcMatrix(6,:),fMatrix(6,:)/gcList(6),'x'); 
% scatter(ntgcMatrix(8,:),fMatrix(8,:)/gcList(8),'p'); 
% hold off;
% legend({'\Gamma_C = 0.02'
%         '\Gamma_C = 0.04'
%         '\Gamma_C = 0.06'
%         '\Gamma_C = 0.08'});
% xlabel('\Phi\tau^{2}\Gamma_c','FontSize', 20);
% ylabel('\Delta \nu/\Gamma_c','FontSize', 20);
% ax = gca;
% ax.FontSize = 20;
% ax.YScale = 'log';
% ax.XLim = [0 100];
% ax.YLim = [0.02 inf];
% 
% fprintf('Linewidth shown. Finished.\n');
% pause;
%}
%% g(2)(0) function

g2Matrix = zeros(nMaxGc,nMaxDens);
for j = 1:nMaxDens
    for i = 1:nMaxGc
        %get q and p for a single simulation
        q(:,:) = -JyMatrix(i,j,:,n0_nStore:nStore);
        p(:,:) = JxMatrix(i,j,:,n0_nStore:nStore);
        g2 = 0;
        for k = 1:m
            %g2(0)
            q4 = mean(q(:,k).^4);
            p4 = mean(p(:,k).^4);
            q2p2 = mean(q(:,k).^2.*p(:,k).^2);
            q2 = mean(q(:,k).^2);
            p2 = mean(p(:,k).^2);
            g2 = g2 + (q4+p4+2*q2p2-8*(q2+p2)+8)/(q2+p2-2)^2;
        end
        g2Matrix(i,j) = g2/m;
    end
end

g2List = reshape(g2Matrix, [1, ntgcSize]);
g2List = g2List(ntgcOrder);

%plot data
figure(3);
hold on;
scatter(ntgcList,g2List,'+');
% scatter(ntgcList,f2List,'o');
xlabel('\Phi\tau^{2}\Gamma_c','FontSize', 20);
ylabel('g^{(2)}(0)','FontSize', 20);
ax = gca;
ax.FontSize = 20;
% ax.YScale = 'log';
% ax.XLim = [0 100];
% ax.YLim = [0 3];

fprintf('Program paused. Press enter to continue.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intensity_print = [[0,gcList]', [densList;intensitySS_mean]];
% save vz=0.1_intensity.dat intensity_print -ascii;
% linewidth_print = [[0,gcList]', [densList;fMatrix]];
% save vz=0.1_linewidth.dat linewidth_print -ascii;
% g2_print = [[0,gcList]', [densList;g2Matrix]];

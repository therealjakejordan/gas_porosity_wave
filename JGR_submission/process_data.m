%% Process data for plots

clear all
close all,
clc
%           
addpath('./Colonel/Hamburger_Helper')
addpath('./Colonel/Smooth_Operators')
%%
tol = 1e-2;
load('./testdata_1.mat')

H           = 5e3;
char        = H/deltam;
T           = ti - ti(1);
dt          = ti(end)-ti(end-1);
indexT      = 1:1:length(ti);
indexZ      = find(GridS.xc < (0.995*char) & GridS.xc > 0.005*char);
%%

phiM1   = phiM(:,indexT)-phi_ss ;
maxphi  = max(phiM1(:));
minphi  = min(phiM1(:));
% umphiXM = umphiXM(:,indexT);
rhophiM = rhophiM(:,indexT);
muM     = muM(:,indexT);
DeltaWM = DeltaWM(:,indexT);
% XM      = umphiXM./(1-phiM1);
rhoM    = rhophiM./phiM1;
%% 
for i = 1:length(GridS.xc)
    [fmc, Pmc, fX, fXccon, im, re] = fft_yang( phiM1(i,1:end-1),dt);
    PmcM(i,:)  = Pmc;
    A(i,:)     = fX;
    Astar(i,:) = fXccon;
    Im(i,:)    = im;
    Re(i,:)    = re;
end
%%
[TC,XCt] = meshgrid(T,GridS.xc);
[FM,XCf] = meshgrid(fmc,GridS.xc);

maxPOWER   = max(PmcM); 
mm         = max(maxPOWER);
normPOWER  = maxPOWER/mm;
ff         = find(normPOWER>0.1/4);
cols       = PmcM(:,ff);
[fpcol, zpcol] = gradient(cols,GridS.dx);
invcols = 0.5 * 1./cols;

qv   = invcols.*zpcol;
qm   = mean(qv,1);
ZK     = repmat(1./GridS.xc,1,length(ff));

save('test_processed.mat','GridS','dt','phiM','phi_ss',...
    'phiM1','PmcM','fmc','ff','A','Astar','Im','Re','TC','XCt','FM','XCf',...
    'maxPOWER','normPOWER','qv','qm','ZK','H','deltam','char','T','indexT','indexZ')


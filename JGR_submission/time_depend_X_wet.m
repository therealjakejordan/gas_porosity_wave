clear all,
close all,
clc
addpath('Colonel/Hamburger_Helper')
addpath('Colonel/Smooth_Operators')
%% Set number of grid points
NV = 2e3;
NS = NV-1;
 
flx = 'gd';
 
tol     = 1e6*eps(1);
theta   = 0.51;
purt    = 0.025;

lambda1  = 1/3;
lambda2  = 1/7;
lambda3  = 1/19; 
lambda4  = 1/23;
lambda5  = 1/29;
lambda6  = 1/31;
lambda7  = 1/53;
lambda8  = 1/67; 
lambda9  = 1/83;
lambda10 = 1/101;
%% Parameters
Wc    = 1;
rhom  = 2500;
rhog0 = 200;
rho0  = rhog0/rhom;
g     = 9.81;
Cg    = 685;
mul0  = 1e9;
xstal = 0.2;
xstalp= 0.6;
einst = 2.5;
suspense = (1 - xstal/xstalp)^(-einst*xstalp);
gamma = 25;
mum0  = mul0 * suspense;
mug   = 1e-5;
k0    = 1e-12;
Ppow  = 0.5;
Pg0n  = (Cg^(2*Ppow)*rhom^Ppow);
s     = 4.11e-6;
Sol   = s*Pg0n;
X0    = Sol*rho0^Ppow;
 
% quantities for scaling
W      = 0.02;  
c0     = mug/k0;
deltam = sqrt(4*mum0/3/c0);
alpha  = c0*W/rhom/g; 
beta   = Cg^2/(g*deltam);
 
H     = 5e3;
k1    = 1e-5;
phi0  = 0.1;
Da    = k1 * (deltam/W);
muw   = mum0*10.^(-gamma*X0);
 
tfrac  = 1e-5/4;
dt     = tfrac*(H/deltam);
 
fprintf('Da = %d\n', Da );
fprintf('delta_0 = %d\n', deltam );
 
%% Gridding.
% Grid where momentum conservation equations will be solved.
Grid.V.xmin = 0;
Grid.V.xmax = H/deltam;
Grid.V.Nx   = NV;
GridV       = build_grid(Grid.V);
% Grid where scalar quantities will live.
Grid.S.xmin = GridV.xc(GridV.dof_xmin);
Grid.S.xmax = GridV.xc(GridV.dof_xmax);
Grid.S.Nx   = NS;
GridS       = build_grid(Grid.S);
%% Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create discrete matrices for Divergance, Gradient and accompanying
% identity matrix
[DV,GV,IV]=build_ops(GridV);
[DS,GS,IS]=build_ops(GridS);
  
%% Initial conditions (Michaut 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a   = phi0/rho0/beta *alpha/(alpha+2*phi0*(1-rho0)); %* alpha/(alpha+2*phi0*(1-rho0));
b   = (1 - phi0 + phi0*rho0)/beta;
c   = (2 - 2*rho0+rho0*alpha) * (phi0/alpha/beta/rho0);
  
wg0 = 1  + (phi0*(1-phi0)*(1-rho0));
phi   = phi0 + (1-phi0)*a*GridS.xc;
rhog  = rho0 * exp(-GridS.xc*b/rho0);
Xl    = X0 * exp(-GridS.xc*c/X0);
X0v   = Xl;
wm    = 1 + a*GridS.xc; wmff = wm;
wg    = wg0*(1 + a*GridS.xc);
  
DeltaW = wm-wg;
umphi  = 1-phi;
umphixstal = umphi*(1-xstal);
rhophi = phi.*rhog;
umphixstalX = umphixstal.*Xl; umphixstalX0 = umphixstalX(1);
  
mu        = ((10.^(-gamma*(X0v))));
 
Gamma     = Da*(Xl-Sol*rhog.^Ppow);
 
rhobc    = rhog(1);
rhophibc = rhog(1)*phi(1);
 
%% Boundary conditions for magma momentum
Param.wm.dof_dir        = [GridS.dof_xmin];
Param.wm.dof_f_dir      = [GridS.dof_f_xmin];      % identify faces on Dirichlet bnd
Param.wm.g              = Wc;                     % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.wm.dof_neu        = [];    % identify cells on Neumann bnd
Param.wm.dof_f_neu      = [];    % identify cells on Neumann bnd
Param.wm.qb             = [];
% Build Boundary conditions
[Bwm,Nwm,fn_wm] = build_bnd(Param.wm,GridS,IS);
%% Boundary conditions for gas mass advection
Param.rhophi.dof_dir        = [GridS.dof_xmin];
Param.rhophi.dof_f_dir      = [GridS.dof_f_xmin];      % identify faces on Dirichlet bnd
Param.rhophi.g              = rhog(1)*phi(1);                 % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.rhophi.dof_neu        = [];    % identify cells on Neumann bnd
Param.rhophi.dof_f_neu      = [];    % identify cells on Neumann bnd
Param.rhophi.qb             = [];
Param.rhophi.dof_out        = GridS.dof_xmax;
% Build Boundary conditions
[Brf,Nrf,fn_rf] = build_bnd(Param.rhophi,GridS,IS);
%% Boundary conditions for gas mass advection
Param.umphi.dof_dir        = [GridS.dof_xmin];
Param.umphi.dof_f_dir      = [GridS.dof_f_xmin];   % identify faces on Dirichlet bnd
Param.umphi.g              = umphi(1);               % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.umphi.dof_neu        = [];    % identify cells on Neumann bnd
Param.umphi.dof_f_neu      = [];    % identify cells on Neumann bnd
Param.umphi.qb             = [];
Param.umphi.dof_out        = GridS.dof_xmax;
% Build Boundary conditions
[Buf,Nuf,fn_uf] = build_bnd(Param.umphi,GridS,IS);
%% Boundary conditions for gas exsolution
Param.xs.dof_dir        = [GridS.dof_xmin];
Param.xs.dof_f_dir      = [GridS.dof_f_xmin];     % identify faces on Dirichlet bnd
Param.xs.g              = umphixstalX(1);            % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.xs.dof_neu        = [];  % identify cells on Neumann bnd
Param.xs.dof_f_neu      = [];  % identify cells on Neumann bnd
Param.xs.qb             = [];
Param.xs.dof_out        = GridS.dof_xmax;
% Build Boundary conditions
[Bxs,Nxs,fn_XS] = build_bnd(Param.xs,GridS,IS);
  
Param.dof_out = GridS.dof_xmax;
%% 
i   = 0;
k   = 0;
f   = 0;
ct  = 0;
 
cnumber = H/deltam;
ndphi = 1;
while ct<=3*H/deltam
    i    = i+1;
    ct   = ct + dt;
    phim = phi;
      
    if isnan(ndphi)==1
        break
    end
         
    if i>1
        wm     = interp1(GridV.xc,wm,GridS.xc,'linear');
        wg     = interp1(GridV.xc,wg,GridS.xc,'linear');
        DeltaW = interp1(GridV.xc,DeltaW,GridS.xc,'linear');
        Gamma  = Da*(Xl-Sol*rhog.^Ppow);
    end
     
     
    if abs(ct/(H/deltam) - 1) < 1e6*eps(1)
       phi_ss = phi; 
       umphixstalX_ss = umphixstalX;
       rhophi_ss     = rhophi;
       mu_ss         = mu;
       DeltaW_ss     = DeltaW;
       wm_ss         = wm; 
    end
     
    if ct > H/deltam 
        dt = tfrac*(H/deltam);
        %Boundary conditions for gas exsolution
        Param.xs.dof_dir        = [GridS.dof_xmin];
        Param.xs.dof_f_dir      = [GridS.dof_f_xmin];     % identify faces on Dirichlet bnd
        Param.xs.g              = umphixstalX;                   % set Dirichlet Bnd
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Param.xs.dof_neu        = [];  % identify cells on Neumann bnd
        Param.xs.dof_f_neu      = [];  % identify cells on Neumann bnd
        Param.xs.qb             = [];
        Param.xs.dof_out        = GridS.dof_xmax;
         
         
         
        Param.xs.g = umphixstalX0*(1 + purt*cos(2*pi*(ct)/(lambda1*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda2*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda3*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda4*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda5*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda6*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda7*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda8*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda9*cnumber))...
            + purt*cos(2*pi*(ct)/(lambda10*cnumber)));
         
        % Build Boundary conditions
        [Bxs,Nxs,fn_XS] = build_bnd(Param.xs,GridS,IS);
    end
     
       
    if ct > H/deltam + 0.1* H/deltam
         
        flx = 'mc';
         
    end
      
    % Magma velocity
    wm(1)=1-a*GridS.dx/2;
    wm(2)=1+a*GridS.dx/2;
      
    for j = 2:(GridS.Nx-1)
          
        fi    =( phi(j+1)+phi (j))/2;
        ro    =(rhog(j+1)+rhog(j))/2;
        drhodz=(rhog(j+1) -rhog(j))/GridS.dx;
          
        rhs(j)=1-fi+fi*ro+beta*drhodz;
          
        DeltaW(j)=(DeltaW(j)+(beta*drhodz+ro)*fi/alpha)/2;
        wg(j) = wm(j) - DeltaW(j);
          
        mm=(1-phi(j)^2)/phi(j);
        pp=(1-phi(j+1)^2)/phi(j+1);
        mm=mu(j)*mm/GridS.dx^2;
        pp=mu(j+1)*pp/GridS.dx^2;
        wm(j+1)=(rhs(j)-mm*wm(j-1)+(mm+pp)*wm(j))/pp;
          
        unmphi(j)=1-phi(j);
        rhophi(j)=rhog(j)*phi(j);
    end
      
    wm(end)      = 2*wm(end-1)-1*wm(end-2);
    wg(end)     = 2*wg(end-1)-1*wg(end-2);
    wg(1)       = wg(2);
    DeltaW(end) = 2*DeltaW(end-1)-1*DeltaW(end-2);
    DeltaW(1)   = wm(1)-wg(1);
      
    unmphi(end)=1-phi(end);
    rhophi(end)=phi(end)*rhog(end);
      
    wm     = interp1(GridS.xc,wm,GridV.xc,'linear','extrap');
    wg     = interp1(GridS.xc,wg,GridV.xc,'linear','extrap');
    DeltaW = interp1(GridS.xc,DeltaW,GridV.xc,'linear','extrap');
      
    % Gas update
     
    Awg        = build_adv_op(wg,rhophi,dt,GS,GridS,Param.rhophi,flx);
    fs_rf      = Gamma;
    L_rhog_I   = (IS + theta*dt*(DS*Awg));
    L_rhog_E   = (IS - (1 - theta)*dt*(DS*Awg));
    RHS_rhog   = dt*(fs_rf + fn_rf) + L_rhog_E*(rhophi);
    rhophi     = solve_lbvp(L_rhog_I, RHS_rhog,Brf,Param.rhophi.g,Nrf);
      
    % Magma update
    Awmup     = build_adv_op(wm,umphi,dt,GS,GridS,Param.umphi,flx); 
    fs_uf     = -Gamma;
    L_uf_I    = IS + theta*dt*(DS*Awmup);
    L_uf_E    = IS - (1 - theta)*dt*(DS*Awmup);
    RHS_umphi = dt*(fs_uf+fn_uf) + L_uf_E*(umphi);
    umphi     = solve_lbvp(L_uf_I,RHS_umphi,Buf,Param.umphi.g,Nuf);
      
    % Xsol update
    Awmxs       = build_adv_op(wm,umphixstalX,dt,GS,GridS,Param.xs,flx); 
    fs_XS       = -Gamma;
    L_XS_I      = IS + theta*dt*(DS*Awmxs);
    L_XS_E      = IS - (1 - theta)*dt*(DS*Awmxs);
    RHS_XS      = dt*(fs_XS+fn_XS) + L_XS_E*umphixstalX;
    umphixstalX = solve_lbvp(L_XS_I, RHS_XS,Bxs,Param.xs.g,Nxs);
      
    Xl   = umphixstalX./(umphi*(1-xstal));
    phi  = 1 - umphi;
    mu   = 10.^(-gamma*Xl);
    rhog = rhophi./phi;
      
    phi(1)=2*phi(2)-phi(3);
    rhog(2) = rho0;

      
    ndphi = norm(phi-phim)./GridS.N;
    mass_sum = sum(rhog.*phi+umphi);
      
     
    if mod(i,1e3)==0
         
        fprintf('t = %d\n', ct/(H/deltam) );
        fprintf('ndphi = %d\n', ndphi );
        fprintf('totalmass = %d\n', mass_sum );
        fprintf('step = %d\n', i );
         

    end
     
    if mod(i,(4e1))==0
         
     if  ct >= (2*H/deltam) 
        k                 = k+1;
        phiM(:,k)         = phi;
        umphixstalXM(:,k) = umphixstalX;
        rhophiM(:,k)      = rhophi;
        muM(:,k)          = mu;
        DeltaWM(:,k)      = DeltaW;
        wmM(:,k)          = wm;
        ti(k)             = ct;
     end
     
    end
     
end
%%
save('testfile_1.mat','GridS','dt','phi','phiM',...
    'phi_ss','umphixstalX_ss','rhophi_ss','mu_ss','DeltaW_ss','wm_ss',...
    'umphixstalXM','rhophiM','muM','DeltaWM','wmM','ti',...
    'lambda1','lambda2','lambda3','lambda4',...
    'lambda5','lambda6','lambda7','lambda8',...
    'lambda9','lambda10','k1','deltam','W','Da','purt','gamma','Cg','mul0','xstal','xstalp')
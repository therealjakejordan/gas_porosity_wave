clear all,
close all,
clc
addpath('./Colonel/Hamburger_Helper')
addpath('./Colonel/Smooth_Operators')
theta = 0;
flx = 'gd';
%% Parameters!
Wc    = 1;
rhom  = 2500;
rhog0 = 200;
rho0  = rhog0/rhom;
g     = 9.81;
Cg    = 685;
mul0  = 1e9;
xstal = 0.5;
xstalp= 0.6;
einst = 2.5;
suspense = (1 - xstal/xstalp)^(-einst*xstalp);
gamma = 0;
mum0  = mul0 * suspense;
mug   = 1e-5;
k0    = 1e-12;
Ppow  = 0.5;
Pg0n  = (Cg^(2*Ppow)*rhom^Ppow);
s     = 4.11e-6;
Sol   = s*Pg0n;
X0    = 0;
Rad   = 5;
 
% quantities for scaling
W      = 0.02;  
c0     = mug/k0;
deltam = sqrt(4*mum0/3/c0);
alpha  = c0*W/rhom/g; 
beta   = Cg^2/(g*deltam);
 
H     = 5e3;
k1    = 1e-4;
phi0  = 0.1;
Da    = k1 * (deltam/W);
muw   = mum0*10.^(-gamma*X0);
 
tfrac  = 1e-4;
dt     = tfrac*(H/deltam);
 
fprintf('Da = %d\n', Da );
fprintf('delta_0 = %d\n', Da );

%% Gridding.
% Need to set up staggered grid.  See Tackley pg. 51
% Grid where momentum conservation equations will be solved.
Grid.xmin = 0;
Grid.xmax = H/deltam;
Grid.N    = 5e2;
Grid.Nx   = Grid.N;
Grid      = build_grid(Grid);

%% Operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create discrete matrices for Divergance, Gradient and accompanying
% identity matrix
[D,G,I]=build_ops(Grid);
dof_f_bnd = [Grid.dof_f_xmin; Grid.dof_f_xmax];
Dm(:,dof_f_bnd) = 0;
%% Initial conditions (Michaut 2013)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a   = phi0/rho0/beta *alpha/(alpha+2*phi0*(1-rho0)); %* alpha/(alpha+2*phi0*(1-rho0));
b   = (1 - phi0 + phi0*rho0)/beta;
c   = (2 - 2*rho0+rho0*alpha) * (phi0/alpha/beta/rho0);
  
wg0 = 1  + (phi0*(1-phi0)*(1-rho0));
% Initial conditions are stolen from Chloe's stuff for point of comparison.
% Make initial conditions
phi   = phi0 + 0*(1-phi0)*a*Grid.xc;
rhog  = rho0 * exp(-Grid.xc*b/rho0);
Xl    = X0 * exp(-Grid.xc*c/X0);
X0v   = Xl;
wm    = (1 + a*Grid.xf); 
wmff  = wm;
wg    = 0*wg0*(1 + a*Grid.xf);
  
DeltaW = wm-wg;
umphi  = 1-phi;
umphixstal = umphi*(1-xstal);
rhophi = phi.*rhog;
umphixstalX = umphixstal.*Xl; umphixstalX0 = umphixstalX(1);
  
mu        = ((10.^(-gamma*(X0v))));
 
Gamma     = 0*Da*(Xl-Sol*rhog.^Ppow);
 
rhobc    = rhog(1);
rhophibc = rhog(1)*phi(1);
%% Rx Boundary conditions for scalar fields
% Boundary conditions for gas mass advection
Param.rhophi.dof_dir        = [Grid.dof_xmin];
Param.rhophi.dof_f_dir      = [Grid.dof_f_xmin];      % identify faces on Dirichlet bnd
Param.rhophi.g              = rhog(1)*phi(1);                 % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.rhophi.dof_neu        = [];    % identify cells on Neumann bnd
Param.rhophi.dof_f_neu      = [];    % identify cells on Neumann bnd
Param.rhophi.qb             = [];
Param.rhophi.dof_out        = Grid.dof_xmax;
% Build Boundary conditions
[Brf,Nrf,fn_rf] = build_bnd(Param.rhophi,Grid,I);
%% Boundary conditions for gas mass advection
Param.umphi.dof_dir        = [Grid.dof_xmin];
Param.umphi.dof_f_dir      = [Grid.dof_f_xmin];   % identify faces on Dirichlet bnd
Param.umphi.g              = umphi(1);               % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.umphi.dof_neu        = [];    % identify cells on Neumann bnd
Param.umphi.dof_f_neu      = [];    % identify cells on Neumann bnd
Param.umphi.qb             = [];
Param.umphi.dof_out        = Grid.dof_xmax;
% Build Boundary conditions
[Buf,Nuf,fn_uf] = build_bnd(Param.umphi,Grid,I);
%% Boundary conditions for magma momentum equation
Param.wm.dof_dir        = [Grid.dof_xmin];
Param.wm.dof_f_dir      = [Grid.dof_f_xmin];      % identify faces on Dirichlet bnd
Param.wm.g              = [Wc];                     % set Dirichlet Bnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Param.wm.dof_neu        = [Grid.dof_xmin];    % identify cells on Neumann bnd
Param.wm.dof_f_neu      = [Grid.dof_f_xmin];    % identify cells on Neumann bnd
Param.wm.qb             = [a];
% Build Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions for matrix velocity
If = speye(Grid.Nf);
B_wm = If(Param.wm.dof_f_dir,:);
N_wm = If; N_wm(:,Param.wm.dof_f_dir) = [];
fn_wm = spalloc(Grid.Nf,1,0);
fn_wm(Param.wm.dof_neu) = Param.wm.qb*Grid.A(Param.wm.dof_f_neu)./Grid.V(Param.wm.dof_neu);
%% 
i   = 0;
k   = 0;
f   = 0;
ct  = 0;
%% 
cnumber = H/deltam;
ndphi = 1;
while ct<=cnumber
    i = i+1;
 
    ct = ct+dt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaging 
    compF  = diag((1-phi.^2)./(phi)); 
    phif   = comp_mean(phi,-1,1,Grid); 
    phif   = diag(phif); 
    phif(1) = phif(2);  phif(end) = phif(end-1);
    rhof   = comp_mean(rhog,-1,1,Grid);
    rhof   = diag(rhof); 
    rhof(1) = rhof(2);  rhof(end) = rhof(end-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gas density gradient with  first/last cells estimated.
    drhodz = G*rhog; drhodz(end) = (4/3 * drhodz(end-1) - drhodz(end-2)/3);
    drhodz(1) = (4/3 * drhodz(2) - drhodz(3)/3);
    % Velocity difference with first/last cells estimated.
    DeltaW = (phif./alpha).*(beta*drhodz+rhof);
    DeltaW(end) = (4/3 * DeltaW(end-1) - DeltaW(end-2)/3);
    DeltaW(1) = (4/3 * DeltaW(2) - DeltaW(3)/3);
    
    % Gas velocity
    wg     = wm - DeltaW;
    % Magma Velocity
    f_wm      = (1-phif).*(1-rhof)/alpha + (beta*drhodz+rhof)/alpha ;
    f_wm(end) = 0;
    RHS_wm    = (f_wm + fn_wm);
    L_wm      = G * compF * D; 
    L_wm(Grid.dof_f_xmax,:) = D(1,:); 
    wm          = solve_lbvp(L_wm, RHS_wm,B_wm,Param.wm.g,N_wm);
    wm(wm<Wc)   = Wc + eps;
    % Gas update
    Awg        = build_adv_op(wg,rhophi,dt,G,Grid,Param.rhophi,flx);
    fs_rf      = Gamma;
    L_rhog_I   = (I + theta*dt*(D*Awg));
    L_rhog_E   = (I - (1 - theta)*dt*(D*Awg));
    RHS_rhog   = dt*(fs_rf + fn_rf) + L_rhog_E*(rhophi);
    rhophi     = solve_lbvp(L_rhog_I, RHS_rhog,Brf,Param.rhophi.g,Nrf);
    % Magma update
    Awmup     = build_adv_op(wm,umphi,dt,G,Grid,Param.umphi,flx);
    fs_uf     = -Gamma;
    L_uf_I    = I + theta*dt*(D*Awmup);
    L_uf_E    = I - (1 - theta)*dt*(D*Awmup);
    RHS_umphi = dt*(fs_uf+fn_uf) + L_uf_E*(umphi);
    umphi     = solve_lbvp(L_uf_I,RHS_umphi,Buf,Param.umphi.g,Nuf);

    phi  = 1 - umphi;
    rhog = rhophi./phi;
    
    if mod(i,1e3)==0
    figure(1)
    i
    subplot(1,5,1)
    plot(phi,Grid.xc)
    
    subplot(1,5,2)
    plot(rhog,Grid.xc)
    
    subplot(1,5,3)
    plot(wm,Grid.xf),
    
    subplot(1,5,4)
    plot(rhog,Grid.xc)

    
    
    drawnow
    
    end
end

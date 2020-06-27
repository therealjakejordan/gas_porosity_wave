%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  , , _    _,  , , _   ,  , _  ,   , _, %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% |\/|'|\  / _ |\/|'|\  | ,|'|\ \  / /_, %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% | `| |-\'\_|`| `| |-\ |/\| |-\ \/`'\_  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% '  ` '  ` _| '  ` '  `'  ` '  `'     ` %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%          '                             %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, 
% close all, 
clc
set(groot,'defaultAxesTickLabelinterpreter','latex')
set(groot,'defaultAxesFontSize',18)
addpath('./Colonel/Hamburger_Helper')
addpath('./Colonel/Smooth_Operators')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms
syms  k om al bet phi0 wg0 rho0 lambda
% Enter the matrix
% Vector of unknowns: [\phi \rho_g  \w_g \w_m]

A11 = (om-k);
A12 = 0;
A13 = 0;
A14 = k*(1-phi0);

A21 = 1i*rho0*(k-om + (k*phi0*(1-phi0)*(1-rho0))/al);
A22 = phi0*( (2*phi0*(1-phi0)*(1-phi0+phi0*rho0))/(al*bet)...
    + 1i*( k-om + (k*(1-2*rho0)*phi0*(1-phi0))/al) );
A23 = -phi0*((1 - phi0 + phi0*rho0)/bet - 1i*rho0*k);
A24 = 0;

A31 = (1-rho0)*(1-phi0);
A32 = -phi0*(1i*bet*k + phi0);
A33 = -al;
A34 = al;

A41 = (1-rho0)*(2*phi0 - 1);
A42 = 0;
A43 = al;
A44 = -al*(1+k^2*(1-phi0^2));

one = collect((A11 * A22 * A33 * A44), om);
two = collect((A11 * A22 * A34 * A43), om);
tre = collect((A14 * A22 * A33 * A41), om);
qua = collect((A14 * A22 * A31 * A43), om);
cin = collect((A11 * A23 * A32 * A44), om);
sei = collect((A14 * A21 * A32 * A43), om);
set = collect((A14 * A23 * A32 * A41), om);

DAme = collect((one-two-tre+qua-cin-sei+set),om);
A = [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44];
% aa= [a   b   c   d;   e   f   g   h;   ii  jj  k   l;   m   n   o   p];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DA  = det(A);
DA = collect(DA,om);
W = solve(DAme==0,om);
W = simplify(W);

root1 = matlabFunction(W(1));
root2 = matlabFunction(W(2));
%%
H = 5e3; 
deltam = 44.267276788012865;
scale  = H/deltam;

%Handle for Sol1 & Sol2: @(al,bet,k,phi0,rho0)
alpha = 10;  p_phi0 = 0.1; p_rho0 = 0.1; p_k = logspace(-2,3,50);

cc     = [0 53 107]/255;
L = 2*pi./p_k;    
Beta   = [1e2 1e3];
linestyle = [':','--'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@(al,bet,k,phi0,rho0)
for i = 1:length(Beta)
    
    Cvec1 = real(root1(alpha,Beta(i),p_k,p_phi0,p_rho0).*L./(2*pi));
    Svec1 = imag(root1(alpha,Beta(i),p_k,p_phi0,p_rho0));
    qvec1 = Svec1./Cvec1;
    omvec1 = Cvec1.*p_k;
    
    Cvec2 = real(root2(alpha,Beta(i),p_k,p_phi0,p_rho0).*L./(2*pi));
    Svec2 = imag(root2(alpha,Beta(i),p_k,p_phi0,p_rho0));
    qvec2 = Svec2./Cvec2;
    omvec2 = Cvec2.*p_k;
    
    [x0,y0,~,~] = intersections(1./L,qvec1,[1./min(L), 1./max(L)],[0 0]);
    
    figure(1)
    subplot(1,2,1)
    semilogx(1./L,qvec1,'linewidth',1.5,'color','k','linestyle',linestyle(i));hold on
    semilogx([1./min(L), 1./max(L)],[0 0],'linewidth',1.5,'linestyle','--','color','k')
    semilogx([x0 x0],[min(qvec1) max(qvec1)],'linewidth',1.5,'linestyle','--','color','k')
    xlabel('$f$','Fontsize',24,'Interpreter','latex')
    ylabel('$\mathcal{C}_1$','Fontsize',24,'Interpreter','latex')
    axis square
    axis tight
    xlim([1e-2 1e1])
    text(1.1e-2,2e-4,'growing waves','Fontsize',14,'Interpreter','latex')
    
    subplot(1,2,2)
    semilogx(1./L,qvec2,'linewidth',1.5,'color','k','linestyle',linestyle(i));hold on
    xlabel('$f$','Fontsize',24,'Interpreter','latex')
    ylabel('$\mathcal{C}_2$','Fontsize',24,'Interpreter','latex')
    axis square
    axis tight
    xlim([1e-2 1e1])
end

fig1=figure(1);
fig1.Renderer='Painters';

legend({'$\beta = 1$e2','$\beta = 1$e3'},'Interpreter','latex', 'FontSize', 18,'Location','southwest')






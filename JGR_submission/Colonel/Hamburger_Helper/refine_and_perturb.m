function [X,Y,Xf,Yf,K,Kf,Kc] = refine_and_perturb(refine,Grid)
% This is gross and needs a lot of work.

%K = squeeze(10*exp(smooth3(smooth3(randn(1,Grid.Nx,Grid.Ny),'gaussian',[3 3 3],sd)...
%    ,'gaussian',[3 3 3],sd)))';
%K = K*sd;

K = squeeze(100*exp(smooth3(smooth3(randn(1,Grid.Nx,Grid.Ny)))))';
% K = 10e3*K;
[X,Y] = meshgrid(linspace(0,Grid.Lx,Grid.Nx),linspace(0,Grid.Ly,Grid.Ny));
[Xf,Yf] = meshgrid(linspace(0,Grid.Lx,refine*Grid.Nx),linspace(0,Grid.Ly,refine*Grid.Ny));

[Xcor,Ycor] = meshgrid(linspace(0,Grid.Lx,Grid.Nx/refine),linspace(0,Grid.Ly,Grid.Ny/refine));
log10Kf = interp2(X,Y,log10(K),Xf,Yf,'spline');

Kcor =squeeze(exp(smooth3(smooth3(randn(1,Grid.Nx/refine,Grid.Ny/refine)))))';
log10Kcor = interp2(Xcor,Ycor,log10(Kcor),X,Y,'spline');
Kf = 10.^log10Kf;
Kc = 10.^log10Kcor;

end


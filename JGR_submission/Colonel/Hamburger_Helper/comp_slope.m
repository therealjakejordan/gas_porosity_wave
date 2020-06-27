function [slope] = comp_slope(c,dx,limiter)

% limiter = indicates type of slope limiter:
%           gd = Godunov 
%           lw = Lax-Wendroff
%           bw = Beam Warming
%           fr = Fromm
%           mm = minmod
%           sb = suberbee
%           mc = MC-limiter (monotenized central-difference limiter)
%           vl = van Leer
%           va = van Alba
%
% Output:
% slope = N by 1 vector of limiter slopes
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> c = sin(Grid.xc);
% >> slope = comp_slope(c,Grid.dx,'mc');
Nx = length(c); slope = zeros(Nx,1);

%% Linear slopes
if strcmp(limiter,'lw') || strcmp(limiter,'bw') || strcmp(limiter,'fr')
    for i=1:Nx
        if strcmp(limiter,'fr')
            slope(2:Nx-1) = (c(3:Nx)-c(1:Nx-2))/2/dx;
        elseif strcmp(limiter,'lw')
            slope(1:Nx-1) = (c(2:Nx)-c(1:Nx-1))/dx;
        elseif strcmp(limiter,'bw')
            slope(2:Nx) = (c(2:Nx)-c(1:Nx-1))/dx;
        end
    end
end

%% Non-linear slopes
if ~strcmp(limiter,'lw') || ~strcmp(limiter,'bw') || ~strcmp(limiter,'fr')
    % Basic linear slopes
    slope_lw = zeros(Nx,1); slope_lw(1:Nx-1) = (c(2:Nx)-c(1:Nx-1))/dx;
    slope_bw = zeros(Nx,1); slope_bw(2:Nx)   = (c(2:Nx)-c(1:Nx-1))/dx;
    slope_fr = zeros(Nx,1); slope_fr(2:Nx-1) = (c(3:Nx)-c(1:Nx-2))/2/dx;
    
    for i=1:Nx
        if strcmp(limiter,'mm')
            slope = min(slope_lw,slope_bw);
        elseif strcmp(limiter,'mc')
            slope = min(slope_fr, min(2*slope_lw, 2*slope_bw));
        elseif strcmp(limiter,'sb')
            slope = max( min(2*slope_lw,slope_bw) , min(slope_lw,2*slope_bw) );
        elseif strcmp(limiter,'vl')
            error('Slope reconstruction for the van Leer limiter is not implemented.')
        elseif strcmp(limiter,'va')
            error('Slope reconstruction for the van Alba limiter is not implemented.')
        end
    end
end
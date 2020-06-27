function [Gcent] = gradient_centered(cent,G,Grid)
% Adjusted gradient for RHS
Gcent     = G*cent;
Gcent(1)  = 2*Gcent(2) - Gcent(3);
Gcent(end)= 2*Gcent(end-1) - Gcent(end-2);
Gcent = interp1(Grid.xf,Gcent,Grid.xc,'linear');
end


function [ rv ] = infinity_bc_incompressible( phi0, beta )
% % Description:
%   For the incompressible case of the BRS
%   sytle compaction equations phi(infinity) must satisfy: 
%   beta*(1-phi) = (phi-phi0)/((phi^2-phi^3)(1-phi0))
%   where beta = deltarho g / c w0(1-phi0)
%   
%   This function takes parameters phi0 and beta and finds the resultant
%   roots of the quartic polynomial obtained by simplifying the above
%   statement.  We are only interested in real answers and phi may not
%   exceed 1 for obvious physical reasons.
%
%   NOTE TO SELF: This root finding method may not be general enough moving
%   forward.  It simply seems to work for obvious values in this simples
%   case
% % Sample Call: only chosen parameter beta and phi0 must be provided,
% %  >> beta = 1e2;
% %  >> phi0 = =0.2
% %  >> rv = infinity_bc_incompressible(phi0, beta)

% intersect root: (- phi0 + z + beta*z^2 - 2*beta*z^3 + beta*z^4, z)

p = [beta -2*beta beta 1 -phi0];
rv = roots(p);
for i = 1:length(rv)
    possible_roots(i) = isreal(rv(i))==1;
end
rv = rv(possible_roots); 
rv(rv>1)  = []; 
rv(rv<0) = [];
end


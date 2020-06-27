function [Gface] = gradient_faces(face,G,Grid)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dface = G*face; dface(1) = 2*dface(2)-dface(3);
dface(Grid.N+1) = 2*dface(Grid.N)-dface(Grid.N-1);
Gface = dface;
end


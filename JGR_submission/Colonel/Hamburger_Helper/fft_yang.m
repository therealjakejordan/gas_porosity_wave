function [f,Pxx, X, Xccon, im, re]=fft_yang(x,dt)
T=1/dt;
X     = fft(x,numel(x))/numel(x);
Xccon = conj(X);
Pxx   = X.* Xccon;
f     = T/numel(x)*(0:round(numel(x)/2));
Pxx   = Pxx(1:numel(f));
im    = imag(X(1:numel(f)));
re    = real(X(1:numel(f)));

% Test if this code is right
% Make an oscillation with 1Hz frequency with w=2*pi and f=2*pi/w=1
    %e.g t=0:0.01:1; dt=0.01 
     %   x=cos(w*t)=cos(2*pi*t); 
     %  run the code, make sure that only non-zero amplitude occurs at f=1.



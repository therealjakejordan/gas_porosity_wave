function c = orangeblue(m)
interp = 100;
% Put in hex code for desired colors
hex = ['#bd5319';'#dddddd';'#00356b'];
vec = [0;interp/2;interp];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = m;
%N = size(get(gcf,'colormap'),1) % size of the current colormap
c = interp1(vec,raw,linspace(interp,0,N),'pchip');
end


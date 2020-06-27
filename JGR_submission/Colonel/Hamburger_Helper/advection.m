function [ x, dt ] = advection(x,v,n,dt,dh,Grid)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
      ddt=dt;
      vp = zeros(Grid.Nx,1);
      vm = zeros(Grid.Nx,1);
      dx = zeros(Grid.Nx,1);
      

      vp(1)=(v(1)+abs(v(1)))/2;
      vm(1)=(v(1)-abs(v(1)))/2;
        
      for i=2:n-1
        vp(i)=(v(i)+abs(v(i)))/2;
        vm(i)=(v(i)-abs(v(i)))/2;
        dx(i)=-(vp(i)*x(i)-vp(i-1)*x(i-1)+vm(i)*x(i+1)-vm(i-1)*x(i))/dh;   
      end

      dx(n)   =2.*dx(n-1) -dx(n-2);
      dx(1)   =2.*dx(2)   -dx(3);

      for i=1:n
      x(i)=x(i)+ddt*dx(i);
      end

      dt=ddt;
end


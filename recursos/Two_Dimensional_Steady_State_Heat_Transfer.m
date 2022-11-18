close all
clear all
clc

n=40;
w=3;
h=3;

x=linspace(0,w,n);
y=linspace(0,h,n);

T=zeros(n);

T(1,1:n)=100 ; %Top
T(n,1:n)=300 ; %Bottom
T(1:n,1)=500 ; %Left
T(1:n,n)=200 ; %Right

tol=1e-6;
error=1;
k=0

while error > tol
    k=k+1
    Told=T
        for i=2:n-1;
            for j=2:n-1;
                T(i,j) = .25 *(T(i,j-1) + T(i-1,j) + T(i+1,j) + T(i,j+1));
            end
        end
        error = max(max(abs(Told-T)));
end

subplot(3,1,1),contour(x,y,T),colormap
title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar
subplot(3,1,2),pcolor(x,y,T),shading interp,
title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar
subplot(3,1,3)
surf(T')
    xlabel('x')
    ylabel('y')
    zlabel('T')
    colorbar
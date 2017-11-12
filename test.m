clc;
close all;
clear all;
%% data
a = 0.93;
b = 5;
t = 2.49;
v = 0.15;
g=2;
h0 =@(x) 1+0.5*(sin(pi*x));
S2 = @(x) (pi/2*(v-1).*cos(pi*(x-t)));
S2Prim = @(x) (0.5*(v-1).*sin(pi*(x-t)));
num = integral(S2,a,b,'AbsTol', 1e-14);
exa = (S2Prim(b)-S2Prim(a));
h = h0(0.25);
u = 0.95;
A = [0 1 ; -u^2+g*h, 2*u];

%% numerical quantites
[S,L]=eig(A);
lambdaMaxNum = max(L(1,1), L(2,2));
lambdaMinNum = min(L(1,1), L(2,2));

%% exact quantities
lambdaMaxExa = u+(g*h)^0.5;
lambdaMinExa = u-(g*h)^0.5;
vMax = 1/(1+(u+(g*h)^0.5)^2)^0.5*[1; u+(g*h)^0.5];
vMin = 1/(1+(-u+(g*h)^0.5)^2)^0.5*[1;u-(g*h)^0.5];
SExa = [vMin, vMax];
LExa = [lambdaMinExa, 0; 0, lambdaMaxExa];

%% check
errMax = A*vMax-lambdaMaxExa*vMax;
errMin = A*vMin-lambdaMinExa*vMin;
ARec = SExa\LExa*SExa;  % with "/" to solve assocated linear system
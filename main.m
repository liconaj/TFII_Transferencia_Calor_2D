clc
close
clear

addpath heattransf2d

T_env = 20;
k_m = 50;
k_g = 0.1;
q = 195e3;
h_env = 10;
cellsize = 4e-2;
celldivisions = 10;
modelo = "modelos/sketch1.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);
heatsystem = heattransf2d(NodeMesh, T_env, [k_m k_g], h_env, q);
heatsystem.showimtemps()

%{
T = [1700 400];    % K
h = [1000 200];    % W/m2K
k = 25;            % W/mK
cellsize = 1e-3;   %m
celldivisions = 10;

modelo = "modelos/bookex.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);
heatsystem = heattransf2d(NodeMesh, T, k, h);
heatsystem.showimtemps()
%}
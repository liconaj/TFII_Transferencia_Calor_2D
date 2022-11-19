clc
close
clear

addpath heattransf2d


k_m = 50;
k_g = 0.1;
q_g = 195e3;
T_env = 20;
h_env = 500;
cellsize = 4e-2;
celldivisions = 10;

modelo = "modelos/sketch1.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);

nodeparams = dictionary();
nodeparams = heattransf2d.setupnk(nodeparams, 0x6, k_m);         % Metal
nodeparams = heattransf2d.setupnk(nodeparams, 0x9, k_g, q_g);    % Camara combustión
nodeparams = heattransf2d.setupnh(nodeparams, 0x1, h_env, T_env);% Ambiente
nodeparams = heattransf2d.setupni(nodeparams, 0x0);              % Aislante

heatsystem = heattransf2d(NodeMesh,nodeparams);
heatsystem = heatsystem.solvesystem();
heatsystem.showimtemps()


%{
cellsize = 1e-3;
celldivisions = 20;
modelo = "modelos/bookex.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);

nodeparams = dictionary();
nodeparams = heattransf2d.setupnk(nodeparams, 0x6, 25);
nodeparams = heattransf2d.setupnh(nodeparams, 0x1, 1000, 1700);
nodeparams = heattransf2d.setupnh(nodeparams, 0xD, 200, 400);
nodeparams = heattransf2d.setupni(nodeparams, 0x0); 

heatsystem = heattransf2d(NodeMesh, nodeparams);
heatsystem = heatsystem.solvesystem();
heatsystem.showimtemps()
%}
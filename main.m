clc
close
clear

addpath heattransf2d

h_env = 10;
T_env = 20;
k = 50;
h1 = 100;
h2 = 100;
hf = 100;
q = 1e3;
A = (0.12*2+0.14) * 0.5; 

cellsize = 2e-2;
celldivisions = 4;
modelo = "modelos/sketch.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);
heatsystem = heattransf2d(NodeMesh);

nodeparams = dictionary();
nodeparams = heattransf2d.setupnk(nodeparams, 0x6, k); % metal
nodeparams = heattransf2d.setupnh(nodeparams, 0x1, h_env, T_env); % ambiente
nodeparams = heattransf2d.setupnh(nodeparams, 0x2, hf, T_env); % aleta
nodeparams = heattransf2d.setupnh(nodeparams, 0xC, h1, T_env); % refrigeracion izq.
nodeparams = heattransf2d.setupnh(nodeparams, 0xD, h2, T_env); % refrigeracion sup.
nodeparams = heattransf2d.setupnq(nodeparams, 0x9, q, A); % flujo calor
nodeparams = heattransf2d.setupni(nodeparams, 0x0); % aislante

heatsystem.setnodeparams(heatsystem, nodeparams);
heatsystem = heatsystem.solvesystem();
heatsystem.showimtemps()
%%

k_m = 50;
k_g = 0.1;
q_g = 195e3;
T_env = 20;
h_env = 10;
cellsize = 4e-2;
celldivisions = 4;

modelo = "modelos/sketch1.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);

nodeparams = dictionary();
nodeparams = heattransf2d.setupnk(nodeparams, 0x6, k_m);         % Metal
nodeparams = heattransf2d.setupnk(nodeparams, 0x9, k_g, q_g);    % Camara combustión
nodeparams = heattransf2d.setupnh(nodeparams, 0x1, h_env, T_env);% Ambiente
nodeparams = heattransf2d.setupni(nodeparams, 0x0);              % Aislante

heatsystem = heattransf2d(NodeMesh);
heatsystem = heatsystem.setnodeparams(nodeparams);
heatsystem = heatsystem.solvesystem();
heatsystem.showimtemps()

%%

cellsize = 1e-3;
celldivisions = 25;
modelo = "modelos/bookex.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);

nodeparams = dictionary();
nodeparams = heattransf2d.setupnk(nodeparams, 0x6, 25);
nodeparams = heattransf2d.setupnh(nodeparams, 0x1, 1000, 1700);
nodeparams = heattransf2d.setupnh(nodeparams, 0xD, 200, 400);
nodeparams = heattransf2d.setupni(nodeparams, 0x0); 

heatsystem = heattransf2d(NodeMesh);
heatsystem = heatsystem.setnodeparams(nodeparams);
heatsystem = heatsystem.solvesystem();
heatsystem.showimtemps()
%%
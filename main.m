clc
close
clear

addpath heattransf2d

u_env = 10; %(m/s) velocidad media del aire
h_env = 2.38 * u_env^0.89;
T_env = 20;
k = 50;
h1 = 50;
h2 = 50;
hf = 100;
q = 1e3;
A = (0.12*2+0.14) * 0.5;

cellsize = 2e-2;
celldivisions = 4;
modelo = "modelos/sketch.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);

heatsystem = heattransf2d(NodeMesh); % Inicialización del sistema

% Especificar tipo y parámetros de cada color del modelo
heatsystem = heatsystem.setupnk(0x6, k); % metal
heatsystem = heatsystem.setupnh(0x1, h_env, T_env); % ambiente
heatsystem = heatsystem.setupnh(0x2, hf, T_env); % aleta
heatsystem = heatsystem.setupnh(0xC, h1, T_env); % refrigeracion izq.
heatsystem = heatsystem.setupnh(0xD, h2, T_env); % refrigeracion sup.
heatsystem = heatsystem.setupnq(0x9, q, A); % flujo calor
heatsystem = heatsystem.setupni(0x0); % aislante

heatsystem = heatsystem.solvesystem(); % resolver el sistem
heatsystem.showimtemps() % mostrar temperaturas

clc
close
clear

addpath heattransf2d

% CONSTANTES
u_env = 10; %(m/s) velocidad media del aire
h_env = 2.38 * u_env^0.89; % convección aire
T_env = 20; % temperatura ambiente
k = 50; % coeficiente de conduccion pieza
A = (0.12*2+0.14) * 0.5;  % Área por la que entra calor
Nu = 3.39; % Número de Nusselt laminar en placas
eD_r = 1e-5; % rugosidad relativa tuberia lisa
v_r = 1e-6; % viscosidad cinematica del agua
Cp_r = 4180; % calor especifico del agua J/Kg°C
k_r = 0.6; %W/m°C conductividad térmica del agua a 20°C
p_r = 997; %kg/m2 densidad del agua
% Refigeracion 1 (izq)
id_r1 = 0xC;
a_r1 = 0.02; b_r1 = 0.16;
A_r1 = a_r1 * b_r1;
Dh_r1 = heattransf2d.calcDh(a_r1,b_r1);
% Refrigeracion 2 (sup)
id_r2 = 0xD;
a_r2 = 0.12; b_r2 = 0.04;
A_r2 = a_r2 * b_r2;
Dh_r2 = heattransf2d.calcDh(a_r2,b_r2);
% Aletas
Z_f = 0.5; % profundidad
Y = 0.24; % altura pared

% VALORES VARIABLES
q = 20e3; % calor que entra a la pieza 
C_r1 = 0.35; %(L/s) velocidad media refrigerante 1
C_r2 = 0.3; %(L/s) velocidad media refrigerante 2
n_f = 100; %número de aletas en pared
L_f = 0.09; %m longitud aleta

Y_f = Y/(2*n_f); % grosor aletas
C_r1 = C_r1 * 1e-3; C_r2 = C_r2 * 1e-3;
u_r1 = C_r1/A_r1; u_r2 = C_r2/A_r2;
Re_r1 = heattransf2d.calcRe(u_r1,Dh_r1,v_r);
Re_r2 = heattransf2d.calcRe(u_r2,Dh_r2,v_r);
mf_r1 = p_r * C_r1;
mf_r2 = p_r * C_r2;
h_1 = heattransf2d.calchref(k_r,Nu,Dh_r1,Re_r1,eD_r);
h_2 = heattransf2d.calchref(k_r,Nu,Dh_r2,Re_r2,eD_r);
[h_f, eta_f, efe_f] = heattransf2d.calcfinheq(k,h_env,Y_f,Z_f,L_f);

% SISTEMA
cellsize = 2e-2;
celldivisions = 4;
modelo = "modelos/sketch.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);
heatsystem = heattransf2d(NodeMesh); % Inicialización del sistema
heatsystem = heatsystem.setupnk(0x6, k); % metal
heatsystem = heatsystem.setupnq(0x9, q, A); % flujo calor
heatsystem = heatsystem.setupni(0x0); % aislante

% RESULTADOS
Z = 0.5; % plano en eje Z a observar
heatsystem = heatsystem.setupnh(0x1, h_env, T_env); % ambiente
heatsystem = heatsystem.setupnh(0x2, h_f, T_env); % aleta
heatsystem = heatsystem.setupnh(0xC, h_1, T_env); % refrigeracion 1
heatsystem = heatsystem.setupnh(0xD, h_2, T_env); % refrigeracion 2
heatsystem = heatsystem.solvesystem(); % resolver el sistema
heatsystem = heatsystem.setTprop(id_r1,mf_r1,Cp_r);
heatsystem.showimtemps(Z) % mostrar temperaturas
Tmax_r1 = heatsystem.getTmax(id_r1,Z);
Tmax_r2 = heatsystem.getTmax(id_r2,Z);
fprintf("Yfin: %0.2f cm\n",Y_f*1e2)
fprintf("h R1: %0.2f \n", h_1)
fprintf("h R2: %0.2f \n", h_2)
fprintf("Tmax R1: %0.2f°C\n",Tmax_r1)
fprintf("Tmax R2: %0.2f°C\n",Tmax_r2)
clc
close
clear

addpath heattransf2d

% CONSTANTES
q = 10e3; %(W) debe ser 10, 20 o 30 para el problema
k = 50; % coeficiente de conduccion pieza
T_env = 20; %(°C) temperatura ambiente
L = 0.5; %(m) largo del sistema
A = (0.12*2+0.14) * L;  % área por la que entra calor

% Ambiente
air = struct();
air.id = 0x1;
air.u = 10; %(m/s) velocidad media
air.p = 1.2041; %(kg/m3) densidad del aire
air.v = 1.8e-5 / air.p; %(m2/s) viscosidad cinematica
air.k = 0.02; %(W/m°C) conductividad térmica a 20-25°C
air.L = L; %(m) longitud pared en la que actúa el aire
air.T = T_env; % temperatura
air.cp = 1000; %(J/kg°C) calor específico del aire

% Agua
wat = struct();
wat.p = 997; %(kg/m3) densidad del agua
wat.v = 1e-6; %(m2/s) viscosidad cinemática
wat.cp = 4186; %(J/kg°C) calor especifico del agua
wat.k = 0.6; %(W/m°C) conductividad térmica del agua a 20°C
wat.T = T_env; %(°C) temperatura de entrada

% Refigeracion 1 (izq)
rf1 = wat;
rf1.id = 0xC; %id color en skecth;
rf1.a = 0.02; %(m) altura canal
rf1.b = 0.16; %(m) base canal
rf1.A = rf1.a * rf1.b; %(m2) área canal
rf1.Nulam = 6.49; % numero de Nusselt laminar con q' constante.

% Refrigeracion 2 (sup)
rf2 = wat;
rf2.id = 0xD; %id color en skecth;
rf2.a = 0.12; %(m) altura canal
rf2.b = 0.04; %(m) base canal
rf2.A = rf2.a * rf2.b; %(m2) área canal
rf2.Nulam = 3.61; % numero de Nusselt laminar con q' constante.

% Aletas
fin = struct();
fin.id = 0x2;
fin.Z = 0.5; %(m) profundidad
fin.Y = 0.24; %(m) altura pared aletas

% VALORES VARIABLES
rf1.C = 0.89e-3; %(m3/s) caudal 
rf2.C = 0.78e-3;
%fin.n = 50; %número de aletas en pared
fin.t = 10e-3; %(m) grosor aletas
fin.L = 10e-2; %(m) longitud de aletas

% Calculos
%fin.t = fin.Y / (2 * fin.n); %(m) grosor aletas
rf1.u = rf1.C / rf1.A; %(m/s) velocidad promedio
rf2.u = rf2.C / rf2.A;
rf1.Pr = heattransf2d.calcPr(rf1.v,rf1.p,rf1.cp,rf1.k); %numero de Prandtl
rf2.Pr = heattransf2d.calcPr(rf2.v,rf2.p,rf2.cp,rf2.k);
air.Pr = heattransf2d.calcPr(air.v,air.p,air.cp,air.k);
rf1.Dh = heattransf2d.calcDh(rf1.a,rf1.b); %(m) diametro hidraulico
rf2.Dh = heattransf2d.calcDh(rf2.a,rf2.b);
rf1.Re = heattransf2d.calcRe(rf1.u,rf1.Dh,rf1.v); %numero de Reynolds
rf2.Re = heattransf2d.calcRe(rf2.u,rf2.Dh,rf2.v);
air.ReL = heattransf2d.calcReL(air.u,air.L,air.v); %numero de Reynolds lineal
rf1.h = heattransf2d.calchint(rf1.k,rf1.Dh,rf1.Pr,rf1.Re,rf1.Nulam); %coeficiente de convección
rf2.h = heattransf2d.calchint(rf2.k,rf2.Dh,rf2.Pr,rf2.Re,rf2.Nulam);
air.h = heattransf2d.calchext(air.k,air.L,air.Pr,air.ReL);
[fin.h, fin.eta, fin.efe] = heattransf2d.calcfinheq(fin.t,fin.L,fin.Z,k,air.h);

% SISTEMA
cellsize = 2e-2;
celldivisions = 6;
modelo = "modelos/sketch.png";
NodeMesh = nodemesh(modelo, cellsize, celldivisions);
heatsystem = heattransf2d(NodeMesh); % Inicialización del sistema
heatsystem = heatsystem.setupnk(0x6, k); % metal
heatsystem = heatsystem.setupnq(0x9, q, A); % flujo calor
heatsystem = heatsystem.setupni(0x0); % aislante
% Configuraciones variables
heatsystem = heatsystem.setupnh(air.id,air.h,air.T); % ambiente
heatsystem = heatsystem.setupnh(fin.id,fin.h,air.T); % aleta
heatsystem = heatsystem.setupnh(rf1.id,rf1.h,rf1.T); % refrigeracion 1
heatsystem = heatsystem.setupnh(rf2.id,rf2.h,rf2.T); % refrigeracion 2
heatsystem = heatsystem.solvesystem(); % resolver el sistema
heatsystem = heatsystem.setTprop(rf1.id, rf1, L);

% RESULTADOS
Z = L; % plano en eje Z a observar
Tprop = heatsystem.getTprop; %aumento temperatura promedio en canales
Tmax_all = heatsystem.getTmax(Z); %temperatura maxima en toda la pieza
Tmax_rf1 = heatsystem.getTmaxc(rf1.id,Z); % temperatura maxima en borde canal 1
Tmax_rf2 = heatsystem.getTmaxc(rf2.id,Z); % temperatura maxima en borde canal 2
Q_rf1 = heatsystem.getHeatConvec(rf1.id); % (W/m2) Flujo de calor en canal 1
Q_rf2 = heatsystem.getHeatConvec(rf2.id); % (W/m2) Flujo de calor en canal 1
Q_fin = heatsystem.getHeatConvec(fin.id); % (W/m2) Flujo de calor en pared aletas
Q_air = heatsystem.getHeatConvec(air.id); % (W/m2)
Q_out = heatsystem.getHeatConvec();
fprintf("\nRESULTADOS\n\n")
fprintf("    Z: %0.1f m\n",Z)
fprintf("Tprop: %0.2f °C/m\n",Tprop)
fprintf(" Q_in: %0.2f W\n",q)
fprintf("Q_out: %0.2f W\n",Q_out)
%fprintf("t_fin: %0.2f cm\n",fin.Y*1e2) %grosor/altura aletas
fprintf("\nCoeficientes de convección\n")
fprintf("h_rf1: %0.2f W/m2°C\n", rf1.h) % coeficiente de convección 1
fprintf("h_rf2: %0.2f W/m2°C\n", rf2.h) % coeficiente de convección 2
fprintf("h_air: %0.2f W/m2°C\n", air.h) % coeficiente de convección 2
fprintf("h_fin: %0.2f W/m2°C\n", fin.h) % coeficiente de conveccion aletas
fprintf("\nFlujo de calor en paredes con convección\n")
fprintf("Q_rf1: %0.2f W\n",Q_rf1)
fprintf("Q_rf2: %0.2f W\n",Q_rf2)
fprintf("Q_air: %0.2f W\n",Q_air)
fprintf("Q_fin: %0.2f W\n",Q_fin)
fprintf("\nVelocidades promedio fluidos\n")
fprintf("u_rf1: %0.2f m/s\n",rf1.u)
fprintf("u_rf2: %0.2f m/s\n",rf2.u)
fprintf("\nTemperaturas máximas\n")
fprintf("    T_max: %0.3f °C\n",Tmax_all)
fprintf("T_max_rf1: %0.2f°C\n",Tmax_rf1)
fprintf("T_max_rf2: %0.2f°C\n",Tmax_rf2)
heatsystem.showimtemps(Z) % mostrar temperaturas
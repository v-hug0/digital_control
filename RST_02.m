%% Atividade 01 - Parte 01
% Descrição da planta
close all
clear
clc
%%
Tsim = 0.005;

% Parâmetros - Buck
Vi = 20;
Vo = 10;
D = Vo/Vi;

f = 10000;
dIL = 0.05;
dVo = 0.1;

L = Vo*(1-D)/(dIL*f);
C = dIL/(8*dVo*f);

Po = 10;
Io = Po/Vo; 
Ro = Vo/Io;
%%

D = 0.5;
% Simulação
sim('justBuck.slx')
% Coleta de dados
t0 = y0.time;
y0 = y0.data;


%%
D = 0.5;
% Resposta ao degrau, diagrama de Bode e plano Z.
Gvd = tf([Vi],[L*C L/Ro 1]);
Ts = 10e-6; 
Hvd = c2d(Gvd,Ts);
%step(Gvv*Vi,'r--',H*Vi,'b-')
step(Gvd)
ylabel('Vo (V)')

figure
bode(Gvd*D)
figure
rlocus(Gvd)
figure
pzplot(Hvd)

%% Atividade 01 - Parte 02

% Identificação de modelo
Tsim = 0.1;
% tr = 5.24e-4;
% tas = 1e-3; % max pulse duration as the settling time
% N = 10; % tamanho do PRBS
% Tb = tas/N;
% Ts = Tb; % tempo de amostragem (s)
fs = 100*1000; % fs = 1000 KHz | fpwm_gen = 10 KHz
Ts = 10e-6; 
%tr = 1.972e-3; % rise time
tr = 1e-3; % rise time
N = 5; % tamanho do PRBS
%p = 21; % submultiple PWMgen frequency
%Tb = p*Ts;
%Tb = 0.000850; % > 197,2us
Tb = 1.3*tr/N;
%N*Tb>tr
%Tb>tr/N
% Simulação
sim('PRBS.slx')
% Coleta de dados
t1 = u.time;
u = u.data;
plot(t1,u,'r') % plot do PRBS (u)
title('PRBS')
hold on
%figure
hold on
t2 = y.time;
y = y.data;
plot(t2,y) % plot da saída (y)
title('Dados de saída')

%% Identificação por mínimos quadrados

% Ordem do modelo
na = 2;
nb = 2;
d = 0;

D = 0.5;

Nsamples = length(y);

phi = [0 0 0 0 0;
       0 0 0 0 0]; % matriz de regressores (phi)
for i = 3:5000
phi(i,:) = [-y(i-1) -y(i-2) u(i) u(i-1) u(i-2)];
end
y_cp = y(1:5000);
theta = inv(phi'*phi)*phi'*y_cp; % cálculo dos parâmetros do modelo
(theta)
a = theta(1:2);
b = theta(3:5);
Gd = tf([b(1) b(2) b(3)],[1 a(1) a(2)],Ts);
figure
%step(Gd*D,'-')
%hold on

D_list = [0.3;0.5;0.7];
D_list = [0.5];
for k=1:length(D_list)
    % Simulação
    D = D_list(k);
    step(Gd*D,'-')
    hold on
    sim('justBuck.slx')
    % Coleta de dados
    t0 = y0.time;
    y0 = y0.data;
    hold on
    plot(t0,y0) % plot da saída (y)
end
ylabel('Vo (V)')
legend('Modelo ARX via MQ','Planta')

% Simulação do modelo adquirido (ymodel)
ymodel = [0;0];
for i = 5000:10001
ymodel(i,1) = [-y(i-1) -y(i-2) u(i) u(i-1) u(i-2)]*theta;
end
figure
plot(t2/Ts,y*D,'b',t2/Ts,ymodel*D,'--r') % plot do modelo adquirido comparado com o modelo simulado
%title('Modelo obtido')
legend('Real','Modelo',Location='best')
% Análise dos resíduos
res = y(5000:10001)-ymodel(5000:10001);
figure
autocorr(res)


%%
Gd = tf([theta(3) theta(4) theta(5)],[1 theta(1) theta(2)],Ts);
[B,A]=tfdata(Gd,'v');
% -------- Regulation Closed loop Design ----------
Tr = 1.8e-3;
MaxOvershoot = 0.0;
% -------- Closed loop performance specification ----------
[Wo,Qsi]=omega_dmp(Tr,MaxOvershoot);
Hc = tf([Wo^2],[1 2*Wo*Qsi Wo^2]);
% -------- Discretization of desired polynomial ----------
Hd = c2d(Hc,Ts,'zoh'); % Discretized desired CL polynomial
[Bd,Ad]=tfdata(Hd,'v');
Polos_d = roots(Ad);
P = poly([Polos_d(1) Polos_d(2)]);
Hs = [1 -1];
Al = conv(A,Hs);
a1 = Al(2);
a2 = Al(3);
a3 = Al(4);
M = [1 0 0 0 ;
a1 1 0 0 ;
a2 0 1 0 ;
a3 0 0 1 ];
P = poly([roots(P)' 0.6]);
x = inv(M)*P';
R = [x(2) x(3) x(4)];
Sl = [1];
S = conv(B(2:3),Sl);
S = conv(S,Hs);
S = S(2:3);
T = P;
% -------- Reference Tracking Specification ----------
[Wo,dmp] = omega_dmp(1.8e-3,0);
H = tf([Wo^2],[1 2*Wo*dmp Wo^2]);
Gd = c2d(H,Ts,'zoh'); % Discretized plant model
[Bm,Am]=tfdata(Gd,'v'); % Polynomials Bm and Am

%%
Ro = Vo/Io;
D = 0.5;
Tsim = 0.03;
% Simulação
sim('justBuckRST_02.slx')
% Coleta de dados

t1_rst1 = u_rst1.time;
u_rst1 = u_rst1.data;
figure
plot(t1_rst1,u_rst1,'b') % plot do sinal de controle
title('Sinal de controle')
t2_rst1 = y_rst1.time;
y_rst1 = y_rst1.data;
figure
plot(t2_rst1,y_rst1,'r') % plot da saída
title('Dados de saída do Buck em malha fechada')


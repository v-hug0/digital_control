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

phi = zeros(na,na);
phi = [0 0 0 0 0;
       0 0 0 0 0]; % matriz de regressores (phi)
for i = 3:5000
phi(i,:) = [-y(i-1) -y(i-2) u(i-1) u(i-2)];
end
y_cp = y(1:5000);
theta = inv(phi'*phi)*phi'*y_cp; % cálculo dos parâmetros do modelo
(theta)
a = theta(1:2);
b = theta(3:5);
Gd = tf([b(1) b(2) b(3) 0],[1 a(1) a(2)],Ts);
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


% Atividade 02 - parte 01 - STANDARD
% Projeto de Controladores Digitais
figure
Ts = 10e-6; 
% -------- Closed loop performance specification ----------
Tr = 1.2e-3;
MaxOvershoot = 0;
[Wo,Qsi]=omega_dmp(Tr,MaxOvershoot);
K=20;
Hc = tf([Wo^2],[1 2*Wo*Qsi Wo^2]);
% ------- Discretization of plant model and desired polynomial
%theta = [-0.6771 -0.0000 0.0000 -0.0000 1.0765];
theta = [-1.873 0.8744 0.0001 0.0152 0.0112];
Gd = tf([theta(3) theta(4) theta(5)],[1 theta(1) theta(2)],Ts);
[B,A]=tfdata(Gd,'v');
d=0;
B = [zeros(1,d) B];
Hd = c2d(Hc,Ts,'zoh'); % Discretized desired CL polynomial
[Bd,Ad]=tfdata(Hd,'v');
Polos_d = roots(Ad);
Pd = poly([Polos_d(1) Polos_d(2)]);
nA=length(A)-1;
nB=length(B)-1;
% Specification of fixed parts Hs and Hs of polynomials R and S,respectively:
Hs=[1 -1]; % Hs - fixed part of S:
Hr=1; % Hr - fixed part of R:
% Building of extended plant polynomials
BB=conv(B,Hr);
AA=conv(A,Hs);
% Definition of the degree of AA and BB
nBB=length(BB)-1-d;
nAA=length(AA)-1;
% Matriz formada pelos elementos do polinômio A:
MA = zeros((nAA+nBB+d),(nBB+d));
for i=1:(nBB+d)
MA(i:nAA+i,i)=AA';
end
% Matriz formada pelos elementos do polinômio B:
MB = zeros((nAA+nBB+d),(nAA));
for i=1:(nAA)
MB(i:nBB+i,i)=BB(d+1:nBB+1+d)';
end
%Formação da matriz M:
MM = [MA MB];
%
%Especificação do pólos auxiliares:
%display(sprintf('Forneça os %d Pólos Auxiliares de malha fechada desejados:', (nAA+nBB+d-3)));
pdf = [-0.001 -0.005];
Pf = [poly(pdf)];
%
%Polinômio característico desejado:
P = conv(Pd,Pf);
%
X = inv(MM)*P';
%
nS = nBB+d-1;
nR = nAA-1;
So=X(1:nS+1);
Ro=X(nS+2:length(X));
R=conv(Hr,Ro)';
S=conv(Hs,So)';
T=sum(R)
Hz=filt(B,A,Ts)
CS=filt(1,S,Ts)
CR=filt(R,1,Ts)
Hrd=CS*Hz
Hmf=T*feedback(Hrd,CR)
config = RespConfig(Amplitude=20);
step(Gd)
hold on
step(Hmf,config)
legend('Malha aberta','Malha fechada')
%ltiview(Hmf)


%%
Ro = Vo/Io;
D = 0.5;
Tsim = 0.03;
% Simulação
sim('justBuckRST_01.slx')
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


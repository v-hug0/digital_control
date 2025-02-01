function [sys,x0,str,ts] = sf_prbs(t,x,u,flag,K,N,Tamost)
% sf_prbs Gera uma seqüência binária pseudo-aleatória.
%
% K - ganho (sinal de saída entre ± K
% N - Número de células (implementado para 5, 7 e 10)
% Tamost - tempo de amostragem

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialização  %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(K,N,Tamost);

  %%%%%%%%%%%%%%%
  % Atualização %
  %%%%%%%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u,K,N);

  %%%%%%%%%%
  % Saídas %
  %%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u,K,N);

  %%%%%%%%%%%%%%%%%%%%%%%%
  % Flags não utilizadas %
  %%%%%%%%%%%%%%%%%%%%%%%%
  case { 1, 4, 9 }
    sys=[];

  %%%%%%%%%%%%%%%%%%%%%%
  % Flags Inexistentes %
  %%%%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

%
%=============================================================================
% mdlInitializeSizes
%=============================================================================

function [sys,x0,str,ts]=mdlInitializeSizes(K,N,Tamost)

global vet
global aux
global pos

vet=ones(1,N);

if N==5,
   pos=[3 5];
elseif N==7,
   pos=[4 7];
elseif N==10,
   pos=[7 10];
else
   error('Número de células inválido!');
   return;
end

aux=xor(vet(1,pos(1,1)),vet(1,pos(1,2)));
vet=[aux vet(:,1:N-1)];

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 0;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [Tamost 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u,K,N)

global vet
global aux
global pos

aux=xor(vet(1,pos(1,1)),vet(1,pos(1,2)));
vet=[aux vet(:,1:N-1)];

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,K,N)

global vet

if vet(1,N)==1,
   sys=[K];
else
   sys=[-K];
end



% end mdlOutputs

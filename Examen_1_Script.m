%% Control Avanzado - Examen: Controladores y Observadores Óptimos

%% Pregunta 1: Solución Ecuación de Riccati

clear; clc;


% Definición de Sistema
A = [22 7 11;-64 -19 -26;-26 -8 -16];
B = [2;-4;-3];
C = eye(length(A));
D = zeros(length(A),length(B(1,:)));

x0 = zeros(3,1);

% Comprobación de Controlabilidad de par (A,B)
Cab = ctrb(A,B);
rCab = rank(Cab); %-> Rango 3 -> Controlable

% Definición de Matrices de Costos Q & R

% Matriz de Costo de Error de Regulación
Q = [158 122 190;122 110 236;190 236 134];

% Matriz de Costo de Acción del Controlador
R = 1;

% Comprobación de Q y R > 0
PDQ = eig(Q);                   %-> No es positiva definida
PDR = eig(R);                   %-> Positiva Definida

% Definición de Matriz de 
syms P11 P22 P33 P12 P13 P23
P = [P11 P12 P13;P12 P22 P23;P13 P23 P33];

% Solución de la Ecuación de Riccati

[PS11,PS12,PS13,PS22,PS23,PS33] = solve((A')*P+P*A+Q-P*B*inv(R)*B'*P == 0,{P11,P12,P13,P22,P23,P33});

% Matrices P de Solución
PS = zeros(3,3,8); d = 0; nd = 0;

for i = 1:1:length(PS11)
    
    % Matrices P soluciones a Ecuación de Riccati

    PS(:,:,i) = [PS11(i) PS12(i) PS13(i);...
                 PS12(i) PS22(i) PS23(i);...
                 PS13(i) PS23(i) PS33(i)];
    
    % Evaluación por Menores Principales 
    MP = PS(:,:,i);
    MP1 = det(MP(1:1));
    MP2 = det(MP(1:2,1:2));
    MP3 = det(MP(1:3,1:3));
    
    % Evaluación por Valores Característicos
    eigPS = real(eig(PS(:,:,i)));
    
    % Evaluación de valores característicos del sistema en lazo cerrado
    KPS = inv(R)*transpose(B)*PS(:,:,i);
    LambdaPS = eig(A-B*KPS);

    if LambdaPS(1) < 0 && LambdaPS(2) < 0 && LambdaPS(3) < 0 %real(eigPS(3))>= 0 && imag(eigPS) == 0 && real(eigPS(1))>= 0 && real(eigPS(2))>= 0
        d = d + 1;
        PSD(:,:,d) = PS(:,:,i);
    else
        nd = nd + 1;
        PND(:,:,nd) = PS(:,:,i);
    end

end

% Número de Soluciones para P

NumSol = max(d);

% Matriz P Solución de Ecuación de Riccati

P1 = PSD;

% Matriz de Ganancias 

K1 = inv(R)*transpose(B)*PSD;

% Valores característicos del Sistema en Lazo Cerrado

Lambda1 = eig(A-B*K1);

% LQR Comando Matlab

[K,S,CLP] = lqr(A,B,Q,R);

K2 = K; P2 = S; Lambda2 = CLP;

% Solución Ecuación de Riccati por Método Numérico de Euler

Pk = zeros(3)*0.01; 
h = 0.01;

for i = 1:1000
    Pk = Pk + h*(A'*Pk + Pk*A + Q - Pk*B*inv(R)*B'*Pk);
end

P3 = Pk; K3 = inv(R)*transpose(B)*P3; Lambda3 = eig(A-B*K3);

%% Pregunta 2: Filtro de Kalman
clear; clc; 

% Sistema Lineal
A = [-5 -4 -4;0 -3 -2;0 2 -3];
B = [-2;0;1];
C = [1 0 4];
D = 0;

% Perturbaciones y Ruido de Medición
E = [-2;0;1];

varxi  = (0.62809)^2;
vareta = (0.50345)^2;

Sxi  = diag(varxi);
Seta = diag(vareta);

% Comprobaciones
o = rank(obsv(A,C));    % Evaluación de Observabilidad par (A,C)
xi = eig(E*Sxi*E');     % Evaluación si SigmaXi es P.D
eta = eig(C'*Seta*C);   % Evaluación si SigmaEta es P.D

% Condición Inicial, Paso y No. de Iteraciones del Método Numérico de Euler

Pk = eye(3)*0.2572; 
h = 0.02848;
itr = 5;

% Solución por Método Númerico de Ecuación de Riccati para Filtros

for i = 1:itr
    Pk = Pk + h*(A*Pk + Pk*A' + E*Sxi*E' - Pk*C'*inv(Seta)*C*Pk);
end

% Matriz de Ganancias del Observador
P1 = Pk;
L1 = Pk*C'*inv(Seta);
E1 = eig(A-L1*C); 

% Solución Ecuación de Riccati para Filtros por Comando LQR de Matlab

L2 = lqr(A',C',E*Sxi*E',Seta)';


% Solución Ecuación de Riccati para Filtros por Comando LQE de Matlab

Q = Sxi;
R = Seta;
G = E;

% dot(x) = Ax + Bu + Gw            {State equation}
%      y = Cx + Du + v             {Measurements}

[L3,P3,E3] = lqe(A,G,C,Q,R);


%% Pregunta 3: Control LQG

clear;
clc;

% Sistema
A = [-13 8 8;-28 11 16;8 0 -5];
B = [-1;-2;0];
C = [6 -4 -5];
D = 0;

sys = ss(A,B,C,D);

% Perturbaciones y Ruido de Medición
E = [-1;0;-1]; G = E;

varxi  = (0.46300)^2; Xi = varxi;
vareta = (0.33238)^2; Eta = vareta;

Qn = E*Xi*E';
Rn = C*Eta*C';

% Costo de Error de Regulación y Energía de Control
Qc = [66 29 86;29 53 89;86 89 179];
Rc = 4;

% Comprobaciones
o = rank(obsv(A,C));    % Evaluación de Observabilidad par (A,C)
c = rank(ctrb(A,B));    % Evaluación de Controlabilidad par (A,B)
xi = eig(E*varxi*E');     % Evaluación si SigmaXi es P.D
eta = eig(C'*vareta*C);   % Evaluación si SigmaEta es P.D

% Control LQG

[K_3,S_3,CLP_3] = lqr(A,B,Qc,Rc);

[L_3,P_3,E_3] = lqe(A,G,C,Xi,Eta);

% Sistema de Retroalimentación por Estados Estimados en Lazo Cerrado

Lambda_3 = eig(A - B*K_3 - L_3*C);

% Comando LQG Matlab

QWV = blkdiag(Qn,Rn);
QXU = blkdiag(Qc,Rc);

KLQG = lqg(sys,QXU,QWV);


%% Pregunta 4: Control LQG de Sistema No Lineal

clear all; close all; clc;

% Tiempo de Valor de Salida Solicitado

t1 = 0.4787;

% Condiciones Iniciales para el Controlador

x0  = [-5;-3;1];
x0b = [0;0;0]; 

% Punto de Equilibrio a Estabilizar

ye = 8;

% Sistema
A = [0 -1 3;-2 -4 -2;-4 3 -7];
B = [-2;1;0];
C = [8,4,5];
D = 0;

% Comprobaciones
o = rank(obsv(A,C));    % Evaluación de Observabilidad par (A,C)
c = rank(ctrb(A,B));    % Evaluación de Controlabilidad par (A,B)

% Cambio de Coordenadas 
Ae = [A B;C D] ;
Be = [0;0;0;ye];

% Sistema en Punto de Equilibrio
Xe = Ae\Be;

xe = Xe(1:3);
ue = Xe(4);


% Perturbaciones y Ruido de Medición
E = [-2;0;3];       G = E;

VarXi  = 0.53157;   Qn = VarXi^2;
VarEta = 0.52173;   Rn = VarEta^2;

% Costo de Error de Regulación y Energía de Control

Qc = [44 18 40;18 27 36;40 36 56];
Rc = 3;

% Cálculo de Matriz de Ganancias de Controlador K y Observador L

[L,P,E]   = lqe(A,G,C,Qn,Rn);

[K,S,CLP] = lqr(A,B,Qc,Rc);


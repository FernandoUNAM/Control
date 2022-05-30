function [P1,P2,P3,K1,K2,K3,L1,L2,L3] = LQR(A,B,Q,R,Pk,h,itr)

%% Solución Ecuación de Riccati por Solve y Funciones Manuales

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

L1 = eig(A-B*K1);

%% Solución Ecuación de Riccati por Comando Matlab

[K,S,CLP] = lqr(A,B,Q,R);

K2 = K; P2 = S; L2 = CLP;

%% Solución Ecuación de Riccati por Método Numérico de Euler

for i = 1:itr
    Pk = Pk + h*(A'*Pk + Pk*A + Q - Pk*B*inv(R)*B'*Pk);
end

P3 = Pk; K3 = inv(R)*transpose(B)*P3; L3 = eig(A-B*K3);


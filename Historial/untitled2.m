%% Control Avanzado - Examen: Controladores y Observadores Óptimos

%% Pregunta 1

clear; clc;

% Definición de Sistema
A = [22 7 11;-64 -19 -26;-26 -8 -16];
B = [2;-4;-3];

% Comprobación de Controlabilidad de par (A,B)
Cab = ctrb(A,B);
rCab = rank(Cab); %-> Rango 3 -> Controlable

% Definición de Matrices de Costos Q & R

% Matriz de Costo de Error de Regulación
Q = [158 122 190;122 110 236;190 236 134];

% Matriz de Costo de Acción del Controlador
R = 1;

% Comprobación de Q y R > 0
det(Q); %-> No es positiva definida
det(R); %-> Positiva Definida

% Definición de Matriz de 
syms P11 P22 P33 P12 P13 P23
P = [P11 P12 P13;P12 P22 P23;P13 P23 P33];

% Solución de la Ecuación de Riccati
PSolve = solve((A')*P+P*A+Q-P*B*inv(R)*B'*P == 0,{P11,P22,P33,P12,P13,P23});

PSolDouble = double(cell2sym(struct2cell(PSolve)));

NoEqns = length(PSolve.P11);
NoVars = length(PSolDouble)/NoEqns;

P_Sol = zeros(6,8);

for i = 1:1:NoEqns
    P_Sol(:,i) = PSolDouble(i:NoEqns:length(PSolDouble),1);
end

%%

for i = 1:1:NoEqns
    for j = 1:1:NoVars
        Current = P_Sol(j,i);
        if j <=3
            Pi(j) = Current;
        else
            Pside(j-3) = Current;
        end
        Pdiag = diag(Pi);
        
    end
end

%%
PSolM = reshape(PSolDouble,8,6);
PS = nan*ones(3,3,NoEqns);

P3 = [PSolM(1,1) PSolM(1,4) PSolM(1,5);...
      PSolM(1,4) PSolM(1,2) PSolM(1,6);...
      PSolM(1,5) PSolM(1,6) PSolM(1,3)]

%%
for i = 1:1:NoEqns
    PS(:,:,i) = [PSolM(i,1) PSolM(i,4) PSolM(i,5);...
                 PSolM(i,4) PSolM(i,2) PSolM(i,6);...
                 PSolM(i,5) PSolM(i,6) PSolM(i,3)]
end

%%
        
        Ps2 = cat(2,zeros(3,1),cat(1,P2,zeros(1,2)));
        Ps3 = toeplitz([0 0 0],[0 0 PSolM(1,5)]);

        Ps = Ps1 + Ps2 + Ps3 + Pi2 + Pi3;
%    end


%%

NoEqns = length(PSolve.P11);
PSol = cell(1,NoEqns); P_PD = cell(1,NoEqns); K = cell(1,NoEqns);
for i = 1:1:NoEqns
    PSol{i} = double([PSolSym(i) PSolSym(i+NoEqns);PSolSym(i+NoEqns) PSolSym(i+(NoEqns*2))]);
    DetP = vpa(det(PSol{i}),3);
    sprintf("El Conjunto de Solución P%d tiene un determinante de %f%+fj",i,real(DetP),imag(DetP))
    double(PSol{i})
    if real(DetP)>= 0 && imag(DetP) == 0
        P_PD{i} = PSol{i};
        if real(DetP)>0 
        sprintf("P%d ES POSITIVA DEFINIDA",i)
        else 
            sprintf("P%d ES POSITIVA SEMIDEFINIDA",i)
        end
    else
        sprintf("P%d NO ES POSITIVA DEFINIDA NI POSITIVA SEMIDEFINIDA",i)
    end
end
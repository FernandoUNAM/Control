%% Ejemplo de Regulación

% Sistema
Am = [-1 0 1;1 -2 0;0 0 -3];
Bm = [0;0;1];
Cm = [1 1 0];

% Salida Deseada
yd = 3;

% Evaluación de Controlabilidad
Cntr = ctrb(Am,Bm);
rCntr = rank(Cntr);

% Cambio de Coordenadas para utilizar Ley de Control 
Ae = [Am Bm;Cm 0];
Be = [0;0;0;yd];

Xe = inv(Ae)*Be;

% Estados de Equilibrio y Entrada de Equilibrio
xe = Xe(1:3);
ue = Xe(4);

% Solution for yd = 3
ke = acker(Am,Bm)

% Vector de Polos Deseados
P = [-5 -6 -7];

% Matriz de Ganancias por medio de Asignación de Polos
k = place(Am,Bm,P);

%Para tener de respuesta los estados
Cstates = eye(3);

% Matriz de Zeros que sea compatible para el Espacio 
Dm = zeros(3,1);

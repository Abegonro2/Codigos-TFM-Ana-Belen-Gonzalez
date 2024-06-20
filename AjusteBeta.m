%% CÁLCULO DEL FACTOR BETA DE UNA EXPONENCIAL A PARTIR DE TRES PUNTOS TEMPORALES
% Ana Belén González Rodríguez-Máster en Física Médica de la Universidad de
% Valencia.
%% Extracción de los datos
%Cargamos un archivo (Datos) con los datos que necesitamos en el siguiente orden:
% Columna 1: Volumen inicial. Columna 2: Tiempo entre la muestra inicial y
% la segunda medida. Columna 3: Volumen de la segunda medida. Columna 4: 
% Tiempo entre la muestra inicial y la tercera medida. Columna 5: Volumen
% en la tercera medida.

%load("datosexcelestudio.mat")
%% Cálculo del factor beta

t1=Datos(:,2);
t2=Datos(:,4);


V0=Datos(:,1);
V1=Datos(:,3);
V2=Datos(:,5);

s=size(V2,1);

syms b;

for i=1:s
    try
    f(b)=(1-(V1(i)/V0(i)).^(-b+1))./(1-(V2(i)./V0(i)).^(-b+1)) -  t1(i)./t2(i);
    betasym(i) = vpasolve(f(b)==0,b);
    catch
        betasym(i)=0;
    end
end
beta=transpose(double(betasym));


%% Comprobación de que la función vpasolve ha devuelto un valor de beta con precisión suficiente 
a=(1-(V1./V0).^(-beta+1))./(1-(V2./V0).^(-beta+1)) -  t1./t2;

if abs(a)>1e-10
    disp('Precisión insuficiente en el cálculo de beta.')
end

%% Cálculo de alpha

for i=1:s
    try
    g(b)=(V1(i)^(-beta(i)+1))/(-beta(i)+1)-(V0(i)^(-beta(i)+1))/(-beta(i)+1)-b*(t1(i));
    alfasym(i) = vpasolve(g(b)==0,b);
    catch
    alfasym(i) =0;
    end
end
alfa=transpose(double(alfasym));
for i=1:s
    if beta(i)==0;
        alfa(i)=nan
    end
end

%% Comprobación de que la función vpasolve ha devuelto un valor de alfa con precisión suficiente


c=(V1.^(-beta+1))./(-beta+1)-(V0.^(-beta+1))./(-beta+1)-alfa.*(t1);

if c>1e-10
    disp('Precisión insuficiente en el cálculo de alfa.')
end
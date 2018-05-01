%% TP 1: Repetidor Analogico vs Repetidor Digital
close all;
clear all;
clc;
%Defino un nuevo colormap para el plot
mycolormap=[lines(7);prism(6);[1 0 1]];
set(groot,'defaultAxesColorOrder',mycolormap);

%% Ejercicio 3
SNR=(-5:1:30);
n=(1:4:25);

%% Probabilidad de error del sistema ANALOGICO
for i=1:length(n)
    for j=1:length(SNR)
        %dB a unidad
        snr=10^(SNR(j)/10);
        
        aux=0;
        for k=0:n(i)-1
            aux=aux+(1+snr^-1)^k;
        end

        snr_rho_n=snr/aux;
        p_ea(i,j)=normpdf(-sqrt(snr_rho_n))/2;
    end
end
%% Probabilidad de error del sistema DIGITAL
for i=1:length(n)
    for j=1:length(SNR)
        p_ed(i,j)=(1-(1-2*qfunc(sqrt(10^(SNR(j)/10))))^n(i))/2;
    end
end

%% Grafico de la probabilidad de error del sistema digital en funcion del SNR
figure;

semilogy(SNR,p_ea(1,:),'DisplayName','P_{e,n=1}^a');
hold on;
grid on;
grid minor;

semilogy(SNR,p_ea(2,:),'DisplayName','P_{e,n=5}^a');
semilogy(SNR,p_ea(3,:),'DisplayName','P_{e,n=9}^a');
semilogy(SNR,p_ea(4,:),'DisplayName','P_{e,n=13}^a');
semilogy(SNR,p_ea(5,:),'DisplayName','P_{e,n=17}^a');
semilogy(SNR,p_ea(6,:),'DisplayName','P_{e,n=21}^a');
semilogy(SNR,p_ea(7,:),'DisplayName','P_{e,n=25}^a');

semilogy(SNR,p_ed(1,:),'DisplayName','P_{e,n=1}^d');
semilogy(SNR,p_ed(2,:),'DisplayName','P_{e,n=5}^d');
semilogy(SNR,p_ed(3,:),'DisplayName','P_{e,n=9}^d');
semilogy(SNR,p_ed(4,:),'DisplayName','P_{e,n=13}^d');
semilogy(SNR,p_ed(5,:),'DisplayName','P_{e,n=17}^d');
semilogy(SNR,p_ed(6,:),'DisplayName','P_{e,n=21}^d');
semilogy(SNR,p_ed(7,:),'DisplayName','P_{e,n=25}^d');


ylim([1e-6 1e0]);
xlim([-5 25]); %Hardcodeado
xlabel('SNR [dB]')
ylabel('P_{e} (SNR) ');
legend('show', 'Location', 'NorthEast');



%% Ejercicio 4, simulación Monte Carlo.
% Comente varias cosas para que lo entiendan ustedes, despues habria que
% revisarlo.
% Fijamos A = 10, h = 0.1 arbitrariamente.

n = 9; % Numero de repetidores
A = 10;
h = 0.1;
res_snr = 1; % Resolucion de la SNR
M = 1e4; % Numero de experimentos por SNR. 10000 son pocos para un grafico hasta 10^-6, pero sino tarda mucho

% Digital

SNR = 5:res_snr:25;
pe_d_sim = zeros (1,length(SNR)); % Inicializacion
for k = 1:length(SNR)
    var = h^2 * A^2 / 10^(SNR(k)/10);
    fallos = 0;
    for i = 1:M
        if (rand(1) > 0.5)
            X = A;
        else X = -A;
        end
        X_dot = X;
        for j = 1:n
            W = randn(1)*sqrt(var); %Normal varianza var.
            if (h*X_dot + W > 0)
                X_dot = A;
            else X_dot = -A;
            end
        end
        if (X ~= X_dot)
            fallos = fallos + 1;
        end
    end
    pe_d_sim(k) = fallos / M; % Promedio. Converge a pe
end

% Analogico

pe_a_sim = zeros (1,length(SNR));
for k = 1:length(SNR)
    G = 1/h * sqrt (10^(SNR(k)/10) / (10^(SNR(k)/10)+1));
    var = h^2 * A^2 / 10^(SNR(k)/10);
    fallos = 0;
    for i = 1:M
        if (rand(1) > 0.5)
            X = A;
        else X = -A;
        end
        X_dot = X;
        for j = 1:n-1
            W = randn(1)*sqrt(var); %Normal varianza var.
            X_dot = (h*X_dot + W) * G;
        end
        W = randn(1)*sqrt(var);
        if (X_dot*h + W > 0)
            X_dot = A;
        else X_dot = -A;
        end
        if (X ~= X_dot)
            fallos = fallos + 1;
        end
    end
    pe_a_sim(k) = fallos / M; % Promedio. Converge a pe
end

%Grafico de comparacion

figure;

semilogy (SNR, pe_d_sim,'DisplayName','P_{e,sim}^d');
hold on;
grid on;
grid minor;
semilogy (SNR, pe_a_sim,'DisplayName','P_{e,sim}^a');

SNR=(-5:1:30);
semilogy (SNR, p_ed(3,:),'DisplayName','P_{e,teo}^d');
semilogy (SNR, p_ea(3,:),'DisplayName','P_{e,teo}^a');

ylim([1e-6 1e0]);
xlim([5 25]);
xlabel('SNR [dB]')
ylabel('P_{e} (SNR) ');
legend('show', 'Location', 'NorthEast');

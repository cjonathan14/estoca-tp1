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
        for k=0:i-1
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
%%


clc
clear all
close all

M = 4; % Número de símbolos en la constelación (por ejemplo, 16-QAM)
s = 100; % Número de símbolos a modular
data = randi([0 M-1],s,1); % Generar secuencia de datos aleatoria

% Definir el tiempo de muestreo y la frecuencia de la portadora
fs = 100000; % Frecuencia de muestreo (Hz)
fc = 1000; % Frecuencia de la portadora (Hz)
l = 4;

qam = qammod(data,M);
% Crear un vector de tiempo para la señal modulada
t_qam = 0 : 1/fs : (length(qam)-1)/fs;

%qam = qammod(data,l');
%t_qam = (0:length(qam)-1)*T;

% Separar partes real e imaginaria de la señal QAM
real_part = real(qam);
imaginary_part = imag(qam);

% Modulación de las partes real e imaginaria con las portadoras
modulated_real = real_part .* cos(2*pi*fc*t_qam);
modulated_imaginary = imaginary_part .* sin(2*pi*fc*t_qam);

% Suma de las señales moduladas para obtener la señal QAM en el dominio del tiempo
qam_time_domain = modulated_real + modulated_imaginary;

demodulated = qamdemod(qam,l)

figure(1)
subplot(3,1,1)
plot(data,'linewidth',2), grid on;
title(' Señal Original');
xlabel("Data");
ylabel("Amplitude");

% Visualizar la señal modulada en el tiempo
figure(1)
subplot(3,1,2)
plot(t_qam, real(qam )); % Parte real de la señal modulada
hold on;
plot(t_qam, imag(qam )); % Parte imaginaria de la señal modulada
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal QAM Modulada en el Tiempo');
legend('Parte Real', 'Parte Imaginaria');
grid on;

figure(1)
subplot(3,1,3)
plot(t_qam,demodulated)
title('Demodulacion QAM');

% Gráfico de la constelación QAM
scatterplot(qam);
title('Constelación QAM');
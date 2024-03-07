
clc;  
clear all; 
close all;  

%% Generar señal 
[x,fs]=audioread('audio2seg.wav'); 
fs = 270000000;
ts=1/fs;
fc = fs/3; 

df=1;
m=x(:,1);                    	
t = (0:ts:((size(m,1)-1)*ts))';  
int_x = cumsum(m)/fs; 

%xfm =cos(2*pi*t(fc+Df*m))=cos(2*pi*fc*t)*cos(2*pi*m*t*Df)-sin(2*pi*fc*t)*sin(2*pi*m*t*Df)
xfm = cos(2*pi*fc*t).*cos(2*pi*fc*int_x*df)-sin(2*pi*fc*t).*sin(2*pi*fc*int_x*df); %FM Modulated signal
% xi=cos(2*pi*fc*int_x); %The i component of carrier wave
% xq=sin(2*pi*fc*int_x) ; %The q component of carrier wave.
t2 = (0:ts:((size(xfm,1)-1)*ts))';  %making time array to do use Hilbert transform .
%t2 = t2(:,ones(1,size(xfm,2)));
%--------Demodulacion con detector de envolvente------------------------%
xfmq = hilbert(xfm).*exp(-j*2*pi*fc*t2);
z = (1/(2*pi*fc))*[zeros(1,size(xfmq,2));diff(unwrap(angle(xfmq)))*fs];
%-----------------------------------
%Plot
%Mensaje original
figure;
subplot(321);  
plot(t,x,'g');
xlabel('Tiempo');
ylabel('Amplitud');
title("Mensaje Original");
hold on;

subplot(323);  
plot(t,xfm);
xlabel('Tiempo ');
ylabel('Amplitud');
title("Senal modulada en FM");
hold on;

subplot(325);
plot(t,z,'b');
xlabel('Tiempo ');
ylabel('Amplitud');
title('Señal demodulada en FM');
ylim([-0.5 0.5]);
hold off;

%----------------------------------------
%Plot frecuencia

% Calcular la Transformada de Fourier de la señal original
N = length(z);
Y = fft(z);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(N/2))/N;

% Gráficos
subplot(3,2,2);
plot(f,P1,'g');
title('Espectro original');
xlabel('Frecuencia (Hz)');

% Calcular la Transformada de Fourier de la señal modulada
N = length(xfm);
Y = fft(xfm);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(N/2))/N;

% Gráficos
subplot(3,2,4);
plot(f,P1);
%plot(f, P1);
title('Espectro modulada');
xlabel('Frecuencia (Hz)');

% Calcular la Transformada de Fourier de la señal modulada
N = length(z);
Y = fft(z);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(N/2))/N;

% Gráficos
subplot(3,2,6);
plot(f,P1,'b');
title('Espectro demodulada');
xlabel('Frecuencia (Hz)');
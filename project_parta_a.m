% Part A: ADAPTIVE NOTCH FILTER.
% define second order IIR notch filter

w=1/3*pi;
w1=(1/2)*pi;
u=0.001; %step size
r=0.9;
N = 3000;
nT =[1:N];
% x=sin(w*nT);
% noise=zeros(1,N);
% noise(N/2)=-1+(1+1)*rand(1,1);
x_org=sin(w*nT)+(0.1+0.1)*rand(1,N);%input sequence+noise
x=x_org+ 5*sin(w1*nT);%sinusoidal interference: 5*sin(w1*nT)

x1=0;
x2=0;
y1=0;
y2=0;
a0=0;
for n=1:N
    e(n)=x(n)+a0*x1+x2; 
    y(n)=e(n)-a0*r*y1-r*r*y2;
    %a(n+1)=a0+u*y(n)*(x1-r*y1);
    %a(n+1)=a0-u*y(n)*x1;
    if abs(a0-u*y(n)*x1) >2 || a0-u*y(n)*x1 == 2
        a(n+1) = 0;
    else
        a(n+1)=a0-u*y(n)*x1;
    end
    a0=a(n+1);
    y2=y1;
    y1=y(n);
    x2=x1;
    x1=x(n);
end

wn = linspace(-pi,pi,N);
den = [1,a(N+1)*r,r*r];
num = [1,a(N+1),1];
[H_ejw,wn] = freqz(num,den,N);
H_ejw_dB = 20*log10(abs(H_ejw));

y=filter(num,den,x);
Y=fft(y,N);

figure('Name','Part(a) Single Channel Adaptive Filter');
subplot(2,3,1)
plot(x(1:100));
title('Input Signal-xn');
xlabel('Samples;n');
ylabel('Magnitude:x');  

subplot(2,3,2)
plot(y(1:100));
title('Output Signal-yn');
xlabel('Samples:n');
ylabel('Magnitude:y');

subplot(2,3,3)
plot(wn,abs(H_ejw));
title('Notch filter');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,4)
x_fft = abs(fft(x,N));
%[x_freq,wx]=freqz(xx,1,512,'whole');
plot(wn,x_fft);
%plot(wx,x_freq);
title('FFT x');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,5)
y_fft = abs(fft(y,N));
%plot(wn,y_fft);
plot(wn,Y);
%[y_freq,wy]=freqz(y,1,512,'whole');
%plot(wy,y_freq);
title('FFT y');
xlabel('Samples');
ylabel('Magnitude');

figure();
nn=1:N;
plot(nn,a(nn));









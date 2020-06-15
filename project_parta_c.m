% Part A: ADAPTIVE NOTCH FILTER.
% define second order IIR notch filter

w=1/3*pi;
w1=(2/3)*pi;
u=0.001; %step size
r=0.9;
N = 3000;
nT =[1:N];
% x=sin(w*nT);
% noise=zeros(1,N);
% noise(N/2)=-1+(1+1)*rand(1,1);
x_org=sin(w*nT)+(0.1+0.1)*rand(1,N);%input sequence+noise
x=x_org+ 5*sin(w1*nT) + 4*sin(0.8*w1*nT);

x1=0;
x2=0;
y1=0;
y2=0;
a0=0;
for n=1:N
    e(n)=x(n)+a0*x1+x2; 
    y(n)=e(n)-a0*r*y1-r*r*y2;
    %a(n+1)=a0+u*y(n)*(x1-r*y1);
    a(n+1)=a0-u*y(n)*x1;
    if abs(a(n+1)) >2 || a(n+1) == 2
        a(n+1) = 0;
    end
    a0=a(n+1);
    y2=y1;
    y1=y(n);
    x2=x1;
    x1=x(n);
end

y1=0;
y2=0;
z1=0;
z2=0;
b0=0;
for n=1:N
    e2(n)=y(n)+b0*y1+y2; 
    z(n)=e2(n)-b0*r*z1-r*r*z2;
    %b(n+1)=b0+u*z(n)*(y1-r*z1);
    b(n+1)=b0-u*z(n)*y1;
    if abs(b(n+1)) >2 || b(n+1) == 2
        b(n+1) = 0;
    end
    b0=b(n+1);
    z2=z1;
    z1=z(n);
    y2=y1;
    y1=y(n);
end

wn1 = linspace(-pi,pi,N);
den1 = [1,a(N+1)*r,r*r];
num1 = [1,a(N+1),1];
[H_ejw1,wn1] = freqz(num1,den1,N);
H_ejw1_dB = 20*log10(abs(H_ejw1));
y=filter(num1,den1,x);
Y=fft(y,N);
wn2 = linspace(-pi,pi,N);
den2 = [1,b(N+1)*r,r*r];
num2 = [1,b(N+1),1];
[H_ejw2,wn2] = freqz(num2,den2,N);
H_ejw2_dB = 20*log10(abs(H_ejw2));
H_ejw = H_ejw1 .* H_ejw2;
z=filter(num2,den2,y);
Z=fft(z,N);

figure('Name','Part(a) Single Channel Adaptive Filter');
subplot(2,3,1)
plot(x);
title('Input Signal');
xlabel('Samples');
ylabel('Magnitude');  

subplot(2,3,2)
plot(z);
title('Output Signal');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,3)
plot(wn,abs(H_ejw));
title('Notch filter');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,4)
x_fft = abs(fft(x,N));
plot(wn,x_fft);
title('FFT x');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,5)
y_fft = abs(fft(y,N));
%plot(wn,y_fft);
plot(wn,Y);
title('FFT y');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,6)
z_fft = abs(fft(z,N));
%plot(wn,y_fft);
plot(wn,Z);
title('FFT z');
xlabel('Samples');
ylabel('Magnitude');


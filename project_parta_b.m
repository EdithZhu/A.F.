%(b)
u=0.001; %step size
r=0.98;
x1=0;
x2=0;
y1=0;
y2=0;
a0=0;
% N = 500;
% w = (pi/3)*([1:N]/N);
% w1 = (2*pi*(1/4 + 0.001.*[1:N]/N)).*([1:N]/N);
% x_org = sin(w)+(0.1+0.1)*randn(1,N); 
% x= x_org + 4*sin(w1);

N = 3000;
nT =[1:N];
w=1/3*pi;
%w1=2*pi*(1/4 + 0.000001*nT);
w1=2*pi*(1/4 + 0.0001*nT);
x_org=sin(w*nT)+(0.1+0.1)*rand(1,N);%input sequence+noise
x=x_org+ 5*sin(w1.*nT);

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


wn = linspace(-pi,pi,N);
den = [1,a(N+1)*r,r*r];
num = [1,a(N+1),1];
[H_ejw,wn] = freqz(num,den,N);
H_ejw_dB = 20*log10(abs(H_ejw));

y=filter(num,den,x);
Y=abs(fft(y,N));

figure('Name','Part(a) Single Channel Adaptive Filter');
subplot(2,3,1)
plot(x);
title('Input Signal');
xlabel('Samples');
ylabel('Magnitude')  

subplot(2,3,2)
plot(y);
title('Output Signal');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,3,3)
plot(wn,H_ejw_dB);
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
%plot(y_fft);
plot(wn, Y);
title('FFT y');
xlabel('Samples');
ylabel('Magnitude');

figure();
nn=1:N;
plot(nn,a(nn));

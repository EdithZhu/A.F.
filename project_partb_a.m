%training sequence s
N=1000;   %sequence length
A=100;    
s=sign(rand(1,N)*2-1);
for m=1:N
    s(m)=s(m)*A;
    if s(m)==-1
        s(m)=-A;
    end
end

%sequence h and s*h
L=4;
h=zeros(1,L+1);
x=zeros(1,N);
%use the given example
h = [0.3, 1, 0.7, 0.3, 0.2];
for m=5:N
    x(m)=h(1)*s(m)+h(2)*s(m-1)+h(3)*s(m-2)+h(4)*s(m-3)+h(5)*s(m-4);
end

%noise sequence n x=s*h+n
n = wgn(1, N, -30, 'dBW');
x=x+n;

%initialize
d=2;
M=20;
hh=zeros(1,M+1);
hh(M/2)=0.5;
y=zeros(1,N);
er=zeros(1,N);
u=0.000001;

%LMS algorithm
for m=100:N
    for i=1:M+1
        y(m)=y(m)+x(m+1-i)*hh(i);
    end
    er(m)=s(m-d)-y(m);
    for i=1:M+1
        hh(i)=hh(i)+2*u*er(m)*x(m+1-i);
    end
end

%plot
figure(1);
plot(s(N-100-d:N-d));
title('Input Signal');
xlabel('Samples');
ylabel('Magnitude')  

hold on;
plot(y(N-100:N));
title('Output Signal');
xlabel('Samples');
ylabel('Magnitude');
hold off;

figure(2);
subplot(2,2,1);
freqz(h,1);   
title('frequency response of the equalizer');

figure(3);
subplot(2,2,2);
freqz(hh,1);
title('frequency response of the channel');

figure(4);
hf=conv(h,hh);
subplot(2,2,3);
freqz(hf,1);
title('frequency response of the equivalent LTI system');


figure(5)
subplot(2,2,1);
zplane(h);
title('pole_zero of the channel');

subplot(2,2,2)
stem(h);
title('impulse response of the channel');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,2,3);
stem(hh);
title('impulse response of the equalizer');
xlabel('Samples');
ylabel('Magnitude');

subplot(2,2,4);
stem(hf);
title('impulse response of the equivalent LTI system');
xlabel('Samples');
ylabel('Magnitude');


figure(6)
plot(er);

















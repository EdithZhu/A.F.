N=100000;
amp=100;
s_rand=sign(rand(1,N)*2-1);
for m=1:N
s_rand(m)=s_rand(m)*amp;
if s_rand(m)== -1
s_rand(m)=-amp;
end
end
% Define sequence h and s_rand*h
L=4;
h=zeros(1,L+1);
x=zeros(1,N);
h(1)=0.3;
h(2)=1;
h(3)=0.7;
h(4)=0.3;
h(5)=0.2;
for m=5:N
x(m)=h(1)*s_rand(m)+h(2)*s_rand(m-1)+h(3)*s_rand(m-2)+h(4)*s_rand(m-3)+h(5)*s_rand(m-4);
end
% Define sequence n (SNR~35dB) and x=s_rand*h+n
n = wgn(1, N, -30, 'dBW');
x=x+n;
% Initialize
M=10;
g=zeros(1,M+1);
g(M/2)=0.1;
y=zeros(1,N);
er=zeros(1,N);
mu=0.00000000001;
% Main loop
for m=100:N
for mm=1:M+1
y(m)=y(m)+x(m+1-mm)*g(mm);
end
er(m)=y(m)*y(m)-amp*amp;
for mm=1:M+1
g(mm)=g(mm)-mu*er(m)*y(m)*x(m+1-mm);
end
end
% Plot result
figure(1);
plot(s_rand(N-100-M/2:N-M/2));
hold
plot(y(N-100:N),'r-');legend('input signal','output signal');
hold off
% figure(2);
% plot(er);
figure(2);
freqz(h,1);
figure(3);
freqz(g,1);
figure;
hf=conv(h,g);
figure(4);
freqz(hf,1);
figure(5);
zplane(h);
figure(6);
stem(h);grid on;xlabel('time');ylabel('Magnitude plot');title('impulse response of h');
figure(7);
stem(g);grid on;xlabel('time');ylabel('Magnitude plot');title('impulse response of the equamplizer upon convergence');
figure(8);
stem(hf);grid on;xlabel('time');ylabel('Magnitude plot');title('impulse response of the equivamplent LTI');
figure(9);
plot(y);
figure(10);%input,output
subplot(1,2,1)
plot(x,'o');title('Input signal');xlabel('n');ylabel('Magnitude')  
subplot(1,2,2)
plot(y,'o');title('Output signal');xlabel('n');ylabel('Magnitude');

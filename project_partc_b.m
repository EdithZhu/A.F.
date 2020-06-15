N=5000;
ww= -pi:2*pi/N:pi;
HH=sqrt(-1)*ww.*(ww<=0.6*pi & ww>=-0.6*pi) + 0.*(ww>=0.6*pi & ww<=pi) +0.*(ww<=-0.6*pi & ww>=-pi);
w= 0:2*pi/N:2*pi;
H=sqrt(-1)*w.*(w<=0.6*pi & w>=0) + 0.*(w>0.6*pi & w<1.4*pi) +sqrt(-1)*(w-2*pi).*(w<=2*pi & w>=1.4*pi);

figure(1);
subplot(211);plot(ww/pi,abs(HH));xlabel('Frequency w'); ylabel('Magnitude plot'); title('Magnitude response of the filter');grid;
subplot(212);plot(ww/pi,angle(HH));xlabel('Frequency w'); ylabel('Phase plot'); title('Phase of the filter');grid;
% subplot(211);plot(w/pi,abs(H));xlabel('Frequency w'); ylabel('Magnitude plot'); title('Magnitude response of the filter');grid;
% subplot(212);plot(w/pi,angle(H));xlabel('Frequency w'); ylabel('Phase plot'); title('Phase of the filter');grid;

h=ifft(H,N);
[HH,W]=DTFT(h,N);
figure(2);
subplot(211);plot(W/pi,abs(HH));xlabel('frequency w'); ylabel('Magnitude plot'); title('frequency response of the filter');grid;
subplot(212);plot(W/pi,angle(HH));xlabel('frequency w'); ylabel('Magnitude plot'); title('frequency response of the filter');grid;


% generate signal 
t=0:N-1;
M  = 25 ;   % filter length M
L=15;
w_seq=0:pi/L:pi;
xn=0;
for w_seq=0:pi/L:pi
    xn = xn + 3*sin(w_seq*t);
end

xnn=xn;%input signal sequence
xn = xn.' ;
dn=conv(h,xnn);
dn=dn(1:N);
dn=dn.';


rho_max = max(eig(xn*xn.'));   % max eigen value of correlation matrix
mu = (1/rho_max) ;    % step size 0 < mu < 1/rho
disp(mu);
[yn,W,en] = lms_algo(xn,dn,M,mu);

% plot input signal
figure;
subplot(3,1,1);
plot(t,xn);grid;ylabel('magnitude plot');xlabel('n');title('input signal');

% plot output signal 
subplot(3,1,2);
plot(t,yn);grid;ylabel('magnitude plot');xlabel('n');title('output signal');

% plot ideal output signal
subplot(3,1,3);
t1=1:N;
plot(t,dn);grid;ylabel('magnitude plot');xlabel('n');title('ideal output signal');



% plot error
figure; 
plot(t,yn,'b',t1,dn,'g',t,dn-yn,'r');grid;
legend('output signal','ideal output signal','error');ylabel('magnitude plot');xlabel('n');title('error analysis');

figure;
X=fft(xnn,N);
Y=fft(dn,N);
fs=2;
f=fs/N*(0:N/2-1);
subplot(221);plot(t,xnn);grid on;xlabel('n');ylabel('magnitude plot');title('orignal signal');
subplot(222);plot(f,abs(X(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of orignal signal');
subplot(223);plot(t,yn);grid on;xlabel('n');ylabel('magnitude plot');title('output signal');
subplot(224);plot(f,abs(Y(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of output signal');


W_impulse=W(:,end).';
[W_freq,WW]=DTFT(W_impulse,512);
figure;
subplot(121);stem(W_impulse);grid on;xlabel('n');ylabel('magnitude plot');title('impulse response of the filter');
subplot(122);plot(WW/pi,abs(W_freq));xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of the filter');
figure;
subplot(211);plot(WW/pi,abs(W_freq));xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of the filter');
subplot(212);plot(WW/pi,angle(W_freq));xlabel('w in pi unit');ylabel('phase plot');grid on;title('freq response of the filter');

%*************************************test case*****************************
x_test = 3*cos(0.1*pi*t)+3*cos(0.2*pi*t)+3*cos(0.3*pi*t);
x_conv_ideal=conv(h,x_test);
x_conv_real=conv(W_impulse,x_test);
figure;
subplot(121);stem(x_conv_ideal(500:700));grid on;xlabel('n');ylabel('magnitude plot');title('impulse response of ideal signal');
subplot(122);stem(x_conv_real(500:700));grid on;xlabel('n');ylabel('magnitude plot');title('impulse response of real signal');
figure;
X_conv_ideal=fft(x_conv_ideal,N);
X_conv_real=fft(x_conv_real,N);
fs=2;
f=fs/N*(0:N/2-1);
subplot(121);plot(f,abs(X_conv_ideal(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of ideal signal');
subplot(122);plot(f,abs(X_conv_real(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of real signal');
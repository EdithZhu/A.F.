% w= -pi:2*pi/1000:pi;
% H=2.*(w<=0.3*pi & w>=-0.3*pi) + 0.*(w<=0.6*pi & w>0.3*pi) + 0.*(w>=-0.6*pi & w<-0.3*pi) + 1.*(w<=pi & w>0.6*pi) + 1.*(w>=-pi & w<-0.6*pi);
% 
% figure();
% plot(w/pi,H);xlabel('Frequency w'); ylabel('Magnitude plot'); title('Magnitude response of the filter');grid;

% generate signal 
N=5000;
t=0:N-1;
M  = 30 ;   % filter length M
L=20;
w_seq=0:pi/L:pi;
xn=1.5;
for w_seq=0:pi/L:pi
    xn = xn + 3*sin(w_seq*t);
end
dn=2*1.5;
for w_seq=0:pi/L:pi
    if w_seq<=0.3*pi
        dn = dn + 2*3*sin(w_seq*t);
    elseif w_seq>=0.6*pi
        dn = dn + 3*sin(w_seq*t);
    end
end
xnn=xn;%input signal sequence
xn = xn.' ;
dn = dn.' ;
        

rho_max = max(eig(xn*xn.'));   % max eigen value of correlation matrix
mu = (1/rho_max) ;    % step size 0 < mu < 1/rho
disp(mu);
[yn,W,en] = lms_algo(xn,dn,M,mu);

% plot input signal
figure;
subplot(3,1,1);
plot(xn(100:300));grid;ylabel('magnitude plot');xlabel('n');title('input signal');

% plot output signal 
subplot(3,1,2);
plot(yn(100:300));grid;ylabel('magnitude plot');xlabel('n');title('output signal');

% plot ideal output signal
subplot(3,1,3);
plot(dn(100:300));grid;ylabel('magnitude plot');xlabel('n');title('ideal output signal');

% plot error
figure; 
plot(t,yn,'b',t,dn,'g',t,dn-yn,'r');grid;
legend('output signal','ideal output signal','error');ylabel('magnitude plot');xlabel('n');title('error analysis');

% plot input output impulse and freq response
figure;
X=fft(xnn,N);
Y=fft(dn,N);
fs=2;
f=fs/N*(0:N/2-1);
subplot(221);plot(t,xnn);grid on;xlabel('n');ylabel('magnitude plot');title('orignal signal');
subplot(222);plot(f,abs(X(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of orignal signal');
subplot(223);plot(t,yn);grid on;xlabel('n');ylabel('magnitude plot');title('output signal');
subplot(224);plot(f,abs(Y(1:N/2)/(N/2)));grid on;xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of output signal');

figure;
W_impulse=W(:,end).';
[W_freq,WW]=DTFT(W_impulse,512);
subplot(311);stem(W_impulse);grid on;xlabel('n');ylabel('magnitude plot');title('impulse response of the filter');
subplot(312);plot(WW/pi,abs(W_freq));xlabel('w in pi unit');ylabel('magnitude plot');grid on;title('freq response of the filter');
subplot(313);plot(WW/pi,angle(W_freq));xlabel('w in pi unit');ylabel('phase plot');grid on;title('freq response of the filter');
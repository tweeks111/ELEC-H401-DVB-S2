clear all; close all;
f = 1;
Fs = 10;
Ts = 1/Fs;
tmin = 0;
tmax = 3;

n1 = tmin/Ts;
n2 = tmax/Ts;
n = n1:n2;
x = cos(2*pi*f*Ts*n);
figure; plot(n*Ts, x);

N1 = 64;
N2 = 128;
N3 = 256;
X1 = abs(fft(x,N1));
X2 = abs(fft(x,N2));
X3 = abs(fft(x,N3));

F1 = [0 : N1 - 1]/N1;
F2 = [0 : N2 - 1]/N2;
F3 = [0 : N3 - 1]/N3;
figure;
subplot(3,1,1)
plot(F1,X1,'-x'),title('N = 64')
subplot(3,1,2)
plot(F2,X2,'-x'),title('N = 128')
subplot(3,1,3)
plot(F3,X3,'-x'),title('N = 256')

X = fftshift(X1);
F = [-N1/2:N1/2-1]/N1;
figure;
plot(F,X1)

% Nb= 1200;                   % Number of bits    
% Nbps= 6;                    % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6)
% CutoffFreq= 1000000;        % CutOff Frequency of the Nyquist Filter
% RollOff= 0.3;               % Roll-Off Factor
% OSF= 2;                     % Oversampling Factor
% Tsymb= 1/(2*CutoffFreq);        % Symbol Period
% BitRate= Nbps/Tsymb;            % Bit Rate
% SymRate= 1/Tsymb;               % Symbol Rate
% Fs = 4*SymRate;                % Sampling Frequency
% 
% N = 501;                                    % N = Number of taps (ODD ONLY)
% df = Fs/N; 
% fmax = df*(N-1)/2;
% dt = 1/Fs;                                 % Delta_t
% 
% i =1;
% for f = -fmax:df:fmax
%    if(abs(f)<=CutoffFreq)
%        H(i)=1;
%    else
%        H(i)=0;
%    end
%    i=i+1;
% end
% 
% fv = -fmax:df:fmax;
% H = ifftshift(H);
% figure;
% plot(fv,H,"-b", fv,H,"b*");
% figure;
% h = ifft(H,"symmetric");
% t = (-(N-1)/2:(N-1)/2)*dt;
% plot(t,h,"-b", t,h,"*b");


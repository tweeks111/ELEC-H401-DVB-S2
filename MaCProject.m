%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc;clear; close all;
%------Parameters------%

Nb= 1200;                   % Number of bits    
Nbps= 6;                    % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6)
CutoffFreq= 1000000;        % CutOff Frequency of the Nyquist Filter
RollOff= 0.3;               % Roll-Off Factor
OSF= 2;                     % Oversampling Factor
T= 1/(2*CutoffFreq);        % Symbol Period
BitRate= Nbps/T;            % Bit Rate
SymRate= 1/T;               % Symbol Rate

%=============================================%
% Mapping
%------------------------

bn_tx = (randi(2,1,Nb)-1)';                 % bn = Binary sequence
In_tx = mapping(bn_tx,Nbps,'qam');          % In = Symbols


fzero = (1+RollOff)/(2*T);                  % fzero = frequency at which the gain is equal to zero
N = 50;                                    % N = Number of samples
fmax = 5*fzero;                             % fmax = max frequency used to calculate the ifft

for i = 1:N/2
    f = i*2*fmax/N;
    if (f<(1-RollOff)/(2*T))
       H(i)=T;
       H(N-i)=T;
    elseif(f<fzero)
       H(i)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2;
       H(N-i)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2;
    else
       H(i)=0;
       H(N-i)=0;
    end
end
H(N)=T;

H=fftshift(H);
h=ifftshift(ifft(H));

fscale = -fmax+2*fmax/N:2*fmax/N:fmax;
tscale = -N/2:1:N/2-1;

figure;
plot(fscale,H,'r-')
hold off;
figure;
plot(tscale,h,'r-');
hold off;

%-------------------------------------%

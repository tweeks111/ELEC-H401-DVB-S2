%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc; clear all; close all;
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
In_tx = mapping(bn_tx,Nbps,'qam');          % In = Symbols sequence at transmitter


N = 10000;                                  % N = Number of samples
k = 10;                                     % k = Number of samples btw every symbol period in h(t)
%fratio = 16;                               % fratio = Ratio between fzero and fmax  -> zero padding
%fzero =(1+RollOff)/(2*T);                  % fzero = frequency at which the gain is equal to zero
dt = T/k;                                   % time scale : 1 sample = dt [s]
df = 1/(dt*N);                              % frequency scale : 1 sample = df [Hz]


for i=1:N
    if ((i-1)*df<(1-RollOff)/(2*T))
       H(i)=T;
    elseif((i-1)*df<(1+RollOff)/(2*T))
       H(i)=T*(1+cos(pi*T*((i-1)*df-(1-RollOff)/(2*T))/RollOff))/2;  
    else
       H(i)=0;
    end
    i = i+1;
end


h = ifftshift(ifft(H,'symmetric'));
tscale = -dt*N/2:dt:(N/2-1)*dt;
fscale = 0: df : N*df-df;

figure;
plot(fscale,H)
hold off;
figure;
plot(tscale,h,'r-');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

es(1:Nb/Nbps)=0;
for t = 1:Nb/Nbps
    for n = 1:Nb/Nbps
        if(t>n)
             es(t)=es(t)+sum(In_tx(n)*h((t-n)*k));
        end 
    end
end

er=es;

y(1:Nb/Nbps)=0;

for t = 1:Nb/Nbps
    for n = 1:Nb/Nbps
        if(t>n)
             y(t)=y(t)+sum(er(n)*h((t-n)*k));
        end 
    end
end

bn_rx = demapping(y',6,'qam');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ber = 0; 
for k=1:length(y)
    if(bn_rx(k) ~= bn_tx(k))
        ber = ber+1;
    end
end
ber = ber/(Nb/Nbps)
%-------------------------------------%

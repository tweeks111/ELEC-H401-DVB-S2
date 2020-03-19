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
In_tx = mapping(bn_tx,Nbps,'qam');          % In = Symbols


N = 5000; %must be even number
i = 1;
fmax = (1+RollOff)/(2*T) + 50000;

for f = 0: fmax/N :fmax
    if (f<(1-RollOff)/(2*T))
       H(i)=T;
    elseif(f<(1+RollOff)/(2*T))
       H(i)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2;  
       H(i)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2; 
    else
       H(i)=0;
       H(i)=0;  
    end
    i = i+1;
end

h = ifft(H);
t = 0:length(h)-1;
figure;
fr = 0: fmax/N :fmax;

plot(fr,H)
hold off;
figure;
plot(t,h,'r-');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

In_rx(1:length(In_tx),1) = 0;

for n=1:length(In_tx)
    for m =1:1:length(In_tx)
        if(n<m)
            In_rx(n,1) = In_rx(n,1) + In_tx(m,1) * h(m-n+1);
        elseif(n>m)
            In_rx(n,1) = In_rx(n,1) + In_tx(m,1) * h(n-m+1);
        else
            In_rx(n,1) = In_rx(n,1) + In_tx(m,1) * h(1);
        end
    end
end 
y = demapping(In_rx,6,'qam');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ber = 0; 
for k=1:length(y)
    if(y(k,1) ~= bn_tx(k,1))
        ber = ber+1;
    end
end
ber = ber/length(bn_tx) %%%%MARCHE PAS...?
%-------------------------------------%

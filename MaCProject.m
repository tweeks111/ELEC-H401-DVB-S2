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


N = 5000;

for i=1:N/2;
    
    f=i*(1+RollOff)/(2*T)/(N/64);

    if (f<(1-RollOff)/(2*T))
       H(N/2+i)=T;
       H(N/2-i+1)=T; 
    elseif(f<(1+RollOff)/(2*T))
       H(N/2+i)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2;  
       H(N/2-i+1)=T*(1+cos(pi*T*(f-(1-RollOff)/(2*T))/RollOff))/2; 
    else
       H(N/2+i)=0;
       H(N/2-i+1)=0;  
    end
end

h=ifftshift(ifft(H));

figure;
plot(1:N,H,'r-')
hold off;
figure;
plot(1:N,h,'r-');
hold off;

%-------------------------------------%

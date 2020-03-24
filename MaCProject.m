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
OSF= 1;                     % Oversampling Factor
Tsymb= OSF*1/(2*CutoffFreq);        % Symbol Period
BitRate= Nbps/Tsymb;            % Bit Rate
SymRate= 1/Tsymb;               % Symbol Rate
Fs = 4*SymRate;                % Sampling Frequency

%=============================================%
% Mapping
%------------------------

bn_tx = (randi(2,1,Nb)-1)';                 % bn = Binary sequence
In_tx = mapping(bn_tx,Nbps,'qam');          % In = Symbols sequence at transmitter


N = 101;                                    % N = Number of taps (ODD ONLY)
df = Fs/N;                                 % Delta_f : 1 tap = df [Hz]
fmax = df*(N-1)/2;
fvector = linspace(-fmax,fmax,N);                   
dt = 1/Fs;                                 % Delta_t
tvector = (-(N-1)/2:(N-1)/2)*dt;

i=1;
for f = -fmax:df:fmax
    
    if (abs(f)<=(1-RollOff)/(2*Tsymb) )
       H(i)=Tsymb;
    elseif(abs(f)<=(1+RollOff)/(2*Tsymb))
       H(i)=Tsymb*(1+cos(pi*Tsymb*(abs(f)-(1-RollOff)/(2*Tsymb))/RollOff))/2;
    else
       H(i)=0;
    end
    i=i+1;
end

f = -fmax:df:fmax;

H = fftshift(H);
h = ifftshift(ifft(H));

 figure;
 plot(f,H, 'r', fvector,H, '*b');
 
hold off;
figure;
plot(tvector, h, 'r', tvector, h, '*b')
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% es(1:Nb/Nbps)=0;
% for t = 1:Nb/Nbps
%     for n = 1:Nb/Nbps
%         if(t>n)
%              es(t)=es(t)+sum(In_tx(n)*h((t-n)*k));
%         end 
%     end
% end
% 
% er=es;
% 
% y(1:Nb/Nbps)=0;
% 
% for t = 1:Nb/Nbps
%     for n = 1:Nb/Nbps
%         if(t>n)
%              y(t)=y(t)+sum(er(n)*h((t-n)*k));
%         end 
%     end
% end
% 
% bn_rx = demapping(y',6,'qam');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ber = 0; 
% for k=1:length(y)
%     if(bn_rx(k) ~= bn_tx(k))
%         ber = ber+1;
%     end
% end
% ber = ber/(Nb/Nbps)
% %-------------------------------------%

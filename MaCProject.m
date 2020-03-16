%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%           John ANTOUN               %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%

%------Parameters------%

Nb=1000;        % Number of bits    
Nbps=6;         % Number of bits per symbol (BPSK=1,QPSK=2,16QAM=4,64QAM=6)
SymRate = 500;       % Symbol rate
BitRate=SymRate*Nbps;      % Bit rate

%=============================================%
% Mapping
%------------------------

bits_tx = randi(2,1,Nb)-1;
symb_tx = mapping(bits_tx,Nbps,"qam");

%-------------------------------------%

%-------------------------------------%
%    Modulation and Coding Project    %
%-------------------------------------%
%   Authors : Theo LEPOUTTE           %
%             John ANTOUN             %
%                                     %
%   Date : March 16, 2020             %
%-------------------------------------%
clc;clear;close all;

blockcodeRate = 0.5;

% 1. Channel encoder  %
%---------------------%

% The P matrix is built from the given H matrix

P = [[0 1 1 1 0];       % P = Parity Check Matrix
     [1 0 1 0 0];
     [1 0 1 0 1];
     [0 0 1 1 1];
     [1 1 0 0 1]].';
 
K = size(P,1);          % Block Vector Length
N = size(P,2)+K;        % Code Vector Length

I = eye(K);             % I = Identity Matrix
H = [I P.']             % H = Parity Check Matrix for systematic code
G = [P I]               % G = Generator Matrix

disp("Orthogonality check :");
disp(mod(G*H',2))       % Orthogonality check 

BER=0;


d = randi(2,1,K)-1      % d = Message Vector

u = mod(d*G,2)          % u = Codeword  -> modulo-2 multiplying
randIndex = randi(N);
r = u;                  % r = Received Vector with an error (r=u+e)

r(randIndex)=~r(randIndex);  
    r
    spaceList=repmat('      ',1,randIndex-1);
    disp(spaceList+"     ^");
    
% 2. Iterative hard decoding  %
%------------------------%

u = hardDecoding(r,H,10)




% 
% errorNb=0;
% for i = 1:N
%      errorNb = errorNb+1;
%    
% end
% errorNb
% if(errorNb>0) BER=BER+1;
% 
% 
% end
% BER=BER/1000


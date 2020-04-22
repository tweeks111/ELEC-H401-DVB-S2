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

P = [[0 1 1 1 0],       % P = Parity Check Matrix
     [1 0 1 0 0],
     [1 0 1 0 1],
     [0 0 1 1 1],
     [1 1 0 0 1]].';
 
K = size(P,1);          % Block Vector Length
N = size(P,2)+K;        % Code Vector Length

I = eye(K);             % I = Identity Matrix
H = [I P.']             % H = Parity Check Matrix for systematic code
G = [P I]               % G = Generator Matrix

disp("Orthogonality check :");
disp(mod(G*H',2))       % Orthogonality check 

d = randi(2,1,K)-1      % d = Message Vector

u = mod(d*G,2)          % u = Codeword  -> modulo-2 multiplying

r = u;                  % r = Received Vector with an error (r=u+e)
    if(r(5)==1) 
        r(5)=0;
    else 
        r(5)=1;
    end
    r
    disp("                             ^");
    
% 2. Iterative hard decoding  %
%------------------------%

hardDecoding(r,H)



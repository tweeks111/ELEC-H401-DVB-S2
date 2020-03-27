function [psd,f] = welch(x,t,L,D)
% INPUTS :
% - x : input signal
% - t : time
% - L : segment size
% - D : delay between segments
% OUTPUTS :
% - psd : power spectral density
% - f : frequency
% Gather segments in a matrix
N = floor((length(x)-L)/D + 1) ;
Xt = [] ;
for nn = 1 :N,
Xt = [Xt , x((nn-1)*D+1 :(nn-1)*D+L).'] ;
end ;
% Multiply with window and compute FFT
W = 1 - ([0 :L-1]-(L-1)/2).^2 * 4/(L+1)^2 ;
Wt = repmat(W',1,N) ;
Yf = fftshift(fft(Xt.*Wt,L,1),1) ;
% Compute PSD
psd = sum(abs(Yf).^2,2)'/N ;
psd = psd/max(psd) ;
f = [-L/2 :L/2-1]/L * 1/(t(2)-t(1)) ;
end

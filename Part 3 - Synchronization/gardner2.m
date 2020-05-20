function [r_corrected,est_time_error] = gardner2(r,K,downsampling_ratio,Tsymb,M)
    
    % NB : the downsampling is done at the same time
    % NB2 : pchip=cubic
    
    y_eps = zeros(1,length(r)/(M/downsampling_ratio));
    est_time_error = zeros(1,length(r)/(M/downsampling_ratio));
    y_eps(1)=r(1);
    
    for n = 1:length(r)/downsampling_ratio-1
        x_vector=((n-1)*M/downsampling_ratio+1:M/downsampling_ratio*n+2);
        t_vector=x_vector*Tsymb*downsampling_ratio/M;
        symbols = r(x_vector);
        t_err = (n*M/downsampling_ratio+1)*Tsymb*downsampling_ratio/M-est_time_error(n);
        t_mid_err = (n*M/downsampling_ratio+1-(M/downsampling_ratio)/2)*Tsymb*downsampling_ratio/M -est_time_error(n);
        y_eps(n+1) = interp1(t_vector,symbols,t_err,'linear');  % => y[n]
        y_mid = interp1(t_vector,symbols,t_mid_err,'linear');   % => y[n-1/2]
        
        est_time_error(n+1)=est_time_error(n)+2*K*real(y_mid*(conj(y_eps(n+1))-conj(y_eps(n))));
    end
    r_corrected = y_eps;
end
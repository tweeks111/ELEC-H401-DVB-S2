function [r_corrected,est_time_error] = gardner(r,K,M)
    
    % NB : the downsampling is done at the same time
    % NB2 : in the course it is written 2*K but in the project explanations
    %       it is written only K. Not really important depending the report
    % NB3 : pchip=cubic
    
    y_eps = zeros(1,length(r)/M);
    est_time_error = zeros(1,length(r)/M);
    y_eps(1)=r(1);
    
    for n = 1:length(r)/M-1
        x_vector=(1+(n-1)*M:M*n+1);
        symbols = r(1+(n-1)*M:M*n+1);
        n_err = n*M+1-est_time_error(n);
        n_mid_err = n*M/2+1 -est_time_error(n);
        y_eps(n+1) = interp1(x_vector,symbols,n_err,'pchip');  % => y[n]
        y_mid = interp1(x_vector,symbols,n_mid_err,'pchip');   % => y[n-1/2]
        
        est_time_error(n+1)=est_time_error(n)+2*K*real(y_mid*(conj(y_eps(n+1))-conj(y_eps(n))));
    end
    r_corrected = y_eps;
end
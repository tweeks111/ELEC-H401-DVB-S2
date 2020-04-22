function [bit_rx] = demapping(symb_rx,Nbps,modulation)

% INPUTS:
% - symb_rx : vector of input symbols (variance 1)
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - bit_rx : vector of ouput bits

Nsymb = size(symb_rx,1); % Number of symbols

switch modulation,
    
    case 'pam'
        
        % Symbol to integer
        sigma = sqrt(sum(([0:2^Nbps-1]-(2^Nbps-1)/2).^2)/2^Nbps); 
        int_rx = sigma * symb_rx + (2^Nbps-1)/2;

        % Integer detection
        int_det = round(int_rx);
        int_det(find(int_det<0)) = 0;
        int_det(find(int_det>2^Nbps-1)) = 2^Nbps-1;

        % Integer to binary
        mapp_rx  = fliplr(de2bi(int_det));

        % Binary to gray
        bit_rx2(:,1) = mapp_rx(:,1);
        for ii = 2:Nbps,
            bit_rx2(:,ii) = xor( mapp_rx(:,ii-1) , mapp_rx(:,ii) );
        end

        bit_rx = reshape(bit_rx2',Nsymb*Nbps,1);

    case 'qam'

        % REAL PART
        NbpsI = Nbps/2; 
        symb_rxI = real(symb_rx);

        % Symbol to integer
        sigmaI = sqrt(sum(([0:2^NbpsI-1]-(2^NbpsI-1)/2).^2)/2^NbpsI); 
        int_rxI = sigmaI * sqrt(2) * symb_rxI + (2^NbpsI-1)/2;

        % Integer detection
        int_detI = round(int_rxI);
        int_detI(find(int_detI<0)) = 0;
        int_detI(find(int_detI>2^NbpsI-1)) = 2^NbpsI-1;

        % Integer to binary
        mapp_rxI  = fliplr(de2bi(int_detI));

        % Binary to gray
        bit_rx2I(:,1) = mapp_rxI(:,1);
        for ii = 2:NbpsI,
            bit_rx2I(:,ii) = xor( mapp_rxI(:,ii-1) , mapp_rxI(:,ii) );
        end

         
        % IMAGINARY PART
        NbpsQ = Nbps/2; 
        symb_rxQ = imag(symb_rx);

        % Symbol to integer
        sigmaQ = sqrt(sum(([0:2^NbpsQ-1]-(2^NbpsQ-1)/2).^2)/2^NbpsQ); 
        int_rxQ = sigmaQ * sqrt(2) * symb_rxQ + (2^NbpsQ-1)/2;

        % Integer detection
        int_detQ = round(int_rxQ);
        int_detQ(find(int_detQ<0)) = 0;
        int_detQ(find(int_detQ>2^NbpsI-1)) = 2^NbpsQ-1;

        % Integer to binary
        mapp_rxQ  = fliplr(de2bi(int_detQ));

        % Binary to gray
        bit_rx2Q(:,1) = mapp_rxQ(:,1);
        for ii = 2:NbpsQ,
            bit_rx2Q(:,ii) = xor( mapp_rxQ(:,ii-1) , mapp_rxQ(:,ii) );
        end
      
         
        % BIT CONCATENATION
        bit_rx = reshape([bit_rx2I,bit_rx2Q]',Nsymb*Nbps,1);
        
end
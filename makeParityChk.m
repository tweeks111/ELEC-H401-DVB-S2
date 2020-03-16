function [c, newH] = makeParityChk(dSource, H, strategy)
% Generate parity check vector bases on LDPC matrix H using sparse LU decomposition
%
%  dSource : Binary source (0/1)
%  H       : LDPC matrix
%  strategy: Strategy for finding the next non-zero diagonal elements
%            {0} First  : First non-zero found by column search
%            {1} Mincol : Minimum number of non-zeros in later columns
%            {2} Minprod: Minimum product of:
%                         - Number of non-zeros its column minus 1
%                         - Number of non-zeros its row minus 1
%           
%  c       : Check bits 
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com

% Get the matric dimension
[M, N] = size(H);
% Set a new matrix F for LU decomposition
F = H;
% LU matrices
L = zeros(M, N - M);
U = zeros(M, N - M);

% Re-order the M x (N - M) submatrix
for i = 1:M

   % strategy {0 = First; 1 = Mincol; 2 = Minprod}
   switch strategy
      
      % Create diagonally structured matrix using 'First' strategy
      case {0}
         
         % Find non-zero elements (1s) for the diagonal
         [r, c] = find(F(:, i:end));
         
         % Find non-zero diagonal element candidates
         rowIndex = find(r == i);
            
         % Find the first non-zero column
         chosenCol = c(rowIndex(1)) + (i - 1);
            
      % Create diagonally structured matrix using 'Mincol' strategy
      case {1}
         
         % Find non-zero elements (1s) for the diagonal
         [r, c] = find(F(:, i:end));
         colWeight = sum(F(:, i:end), 1);
         
         % Find non-zero diagonal element candidates
         rowIndex = find(r == i);
         
         % Find the minimum column weight
         [x, ix] = min(colWeight(c(rowIndex)));
         % Add offset to the chosen row index to match the dimension of the... 
         % original matrix F
         chosenCol = c(rowIndex(ix)) + (i - 1);
             
      % Create diagonally structured matrix using 'Minprod' strategy   
      case {2}
            
         % Find non-zero elements (1s) for the diagonal
         [r, c] = find(F(:, i:end));
         colWeight = sum(F(:, i:end), 1) - 1;
         rowWeight = sum(F(i, :), 2) - 1;
         
         % Find non-zero diagonal element candidates
         rowIndex = find(r == i);
            
         % Find the minimum product
         [x, ix] = min(colWeight(c(rowIndex))*rowWeight);
         % Add offset to the chosen row index to match the dimension of the... 
         % original matrix F
         chosenCol = c(rowIndex(ix)) + (i - 1);
         
      otherwise
         fprintf('Please select columns re-ordering strategy!\n');
      
   end % switch

   % Re-ordering columns of both H and F
   tmp1 = F(:, i);
   tmp2 = H(:, i);
   F(:, i) = F(:, chosenCol);
   H(:, i) = H(:, chosenCol);
   F(:, chosenCol) = tmp1;
   H(:, chosenCol) = tmp2;
                     
   % Fill the LU matrices column by column
   L(i:end, i) = F(i:end, i);
   U(1:i, i) = F(1:i, i);
         
   % There will be no rows operation at the last row
   if i < M           
            
      % Find the later rows with non-zero elements in column i
      [r2, c2] = find(F((i + 1):end, i));          
      % Add current row to the later rows which have a 1 in column i
      F((i + r2), :) = mod(F((i + r2), :) + repmat(F(i, :), length(r2), 1), 2);
                                                           
   end % if
         
end % for i

% Find B.dsource 
z = mod(H(:, (N - M) + 1:end)*dSource, 2);

% Parity check vector found by solving sparse LU
c = mod(U\(L\z), 2); 

% Return the rearrange H 
newH = H;

fprintf('Message encoded.\n');

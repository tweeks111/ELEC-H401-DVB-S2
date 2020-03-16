function H = makeLdpc(M, N, method, noCycle, onePerCol)
% Create R = 1/2 low density parity check matrix
%
%  M        : Number of row
%  N        : Number of column
%  method   : Method for distributing non-zero element
%             {0} Evencol : For each column, place 1s uniformly at random
%             {1} Evenboth: For each column and row, place 1s uniformly at random
%  noCyle   : Length-4 cycle
%             {0} Ignore (do nothing)
%             {1} Eliminate
%  onePerCol: Number of ones per column
%
%  H        : Low density parity check matrix                   
%
%
% Copyright Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com


% Number of ones per row (N/M ratio must be 2)
if N/M ~= 2
   fprintf('Code rate must be 1/2\n');
end
onePerRow = (N/M)*onePerCol;

fprintf('Creating LDPC matrix...\n');

switch method
   % Evencol
   case {0}
      % Distribute 1s uniformly at random within column
      for i = 1:N
         onesInCol(:, i) = randperm(M)';
      end
        
      % Create non zero elements (1s) index
      r = reshape(onesInCol(1:onePerCol, :), N*onePerCol, 1);
      tmp = repmat([1:N], onePerCol, 1);
      c = reshape(tmp, N*onePerCol, 1);
      
      % Create sparse matrix H
      H = full(sparse(r, c, 1, M, N));
      
   % Evenboth
   case {1}
      % Distribute 1s uniformly at random within column
      for i = 1:N
         onesInCol(:, i) = randperm(M)';
      end
        
      % Create non zero elements (1s) index
      r = reshape(onesInCol(1:onePerCol, :), N*onePerCol, 1);
      tmp = repmat([1:N], onePerCol, 1);
      c = reshape(tmp, N*onePerCol, 1);
     
      % Make the number of 1s between rows as uniform as possible     
      
      % Order row index
      [r, ix] = sort(r);
      
      % Order column index based on row index
      for i = 1:N*onePerCol
         cSort(i, :) = c(ix(i));
      end
      
      % Create new row index with uniform weight
      tmp = repmat([1:M], onePerRow, 1);
      r = reshape(tmp, N*onePerCol, 1);
      
      % Create sparse matrix H
      % Remove any duplicate non zero elements index using logical AND
      S = and(sparse(r, cSort, 1, M, N), ones(M, N));
      H = full(S);     
      
end % switch

% Check rows that have no 1 or only have one 1
for i = 1:M
   
   n = randperm(N);
   % Add two 1s if row has no 1
   if length(find(r == i)) == 0
      H(i, n(1)) = 1;
      H(i, n(2)) = 1;
   % Add one 1 if row has only one 1   
   elseif length(find(r == i)) == 1
      H(i, n(1)) = 1;
   end

end % for i

% If desired, eliminate any length-4 cycle
if noCycle == 1
   
   for i = 1:M
      % Look for pair of row - column
      for j = (i + 1):M         
         w = and(H(i, :), H(j, :));
         c1 = find(w);
         lc = length(c1);
         if lc > 1
                       
            % If found, flip one 1 to 0 in the row with less number of 1s
            if length(find(H(i, :))) < length(find(H(j, :)))
               % Repeat the process until only one column left 
               for cc = 1:lc - 1
                  H(j, c1(cc)) = 0;
               end
            else
               for cc = 1:lc - 1
                  H(i, c1(cc)) = 0;
               end
            end % if            
         
         end % if
      end % for j
   end % for i
  
end % if

fprintf('LDPC matrix is created.\n');

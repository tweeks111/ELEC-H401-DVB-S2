

H = [0 1 0 1 1 0 0 1;
     1 1 1 0 0 1 0 0;
     0 0 1 0 0 1 1 1;
     1 0 0 1 1 0 1 0];

bit_piece = [1 1 1 1 0 1 0 1];
result = bit_piece;
[rows,col]=find(H)

[r1,c1]=size(H)

syndrome = mod(result*H.',2)

c_nodes_fs = rows(col==2)

H_nodes = H(c_nodes_fs,:)
[r2,c2]=size(H_nodes);

A = repmat(result,r2,1)

v_nodes_c = H_nodes.*A   %permet d'éliminer la v_node actuelle

nodes_va1 = mod(sum(v_nodes_c,2),2)

diff_val=nodes_va1(nodes_va1>0)
same_val=nodes_va1(nodes_va1==0)

[r3,c3] =  size(diff_val);
[r4,c4] = size(same_val);
B = repmat(result,r3,1);
vote_fn = bitxor(diff_val,B(:,2))
vote_f = [vote_fn;repmat(result(2),r4,1)]
 function v_nodes = hardDecoding(r,H, maxiter)

    % H = N-K*N = r*c
[K,N] = size(H);

v_nodes=r;
v_to_c=zeros(K,N);

syndrome = v_nodes*H';
iter=0;
while(iter<maxiter && norm(mod(syndrome,2))~=0)
    
    %STEP 0 : send v_nodes to c_nodes
    for c = 1:N
       index = find(H(:,c)); 
       v_to_c(index,c)=v_nodes(c);
    end

    %STEP 1 : the message from the c_node to the v_nodes is calculated 
    %          from the previously sent message
    c_to_v=zeros(K,N);
    for k=1:K
        index = find(H(k,:));
        for i = 1:length(index)
            new_index=index;
            new_index(i)=[];
            c_to_v(k,index(i))= mod(sum(v_to_c(k,new_index)),2);
        end
    end
    % STEP 3 : majority voting
    for c=1:N
        index = find(H(:,c));
        vote=(sum(c_to_v(index,c))+v_nodes(c))/(length(index)+1);
        if(vote>0.5)
            v_nodes(c)=1;
        else
            v_nodes(c)=0;
        end   
    end

    syndrome=v_nodes*H.';
    iter = iter +1;
end
end
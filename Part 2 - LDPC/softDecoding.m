function u = softDecoding(r,H,variance,maxiter)

    % H = N-K*N = r*c
[K,N] = size(H);

v_nodes=real(-2*r/variance);
v_to_c=zeros(K,N);

u = ones(1,N);
iter=0;
%STEP 0 : send v_nodes to c_nodes



while(iter<maxiter && norm(mod(u*H',2))~=0)
    
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
            alpha=min(abs(v_to_c(k,new_index)));
            khi=prod(sign(v_to_c(k,new_index)));
            c_to_v(k,index(i))= khi*alpha;
        end
    end
    % STEP 3 : majority voting
    for c=1:N
        index = find(H(:,c));
        vote=v_nodes(c)+sum(c_to_v(index,c));
        for k = 1:length(index)
            v_to_c(index(k),c)=vote-c_to_v(index(k),c);     
        end
        if(vote<0)
            u(c)=1;
        else
            u(c)=0;
        end   
    end
    
    iter = iter +1;
end

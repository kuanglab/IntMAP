function [P,Expression] = EM(RNASeq,PosInfo,maxiter,ThreeSeq,lambda,Trend)
    s = RNASeq - cat(2,0,RNASeq(1:end-1));
    s = s';
    n = length(s);
    q = triu(ones(n,n));
    L = cumsum(PosInfo(:,2)-PosInfo(:,1))';
    
    q = (q./repmat(L,n,1))*max(L);
    P = ones(n,1)/n;
    
    M = diag(ones(1,n)./L);
    if Trend == 1
        E = ThreeSeq(:,2);
    else
        E = flipud(ThreeSeq(:,2));
    end
    for i = 1:maxiter
        P_old = P;
        
        
        alpha = sum((P_old'./L)*sum(s))/sum(E);
        if isinf(alpha)
            alpha = 0;
        end
        
        % E step
        a = (repmat(P',n,1).*q)./(repmat(q*P,1,n)).*repmat(s,1,n);
        
        % M step
        cvx_begin quiet
            cvx_solver sedumi%sdpt3
            cvx_precision('best')
            
            variable x(n)
            minimize( -sum(a*log(x))+ lambda*square_pos(norm(M*x*sum(s)-alpha*E)))
            subject to
            x>=0
            sum(x)==1       
        cvx_end;
        ixx = find(x<=0);
        x(ixx) = 1e-20;
        P = x;
        clear x;
        if max(abs(P-P_old))<1e-5
            break;
        end
    end
    Expression = (P./L')*sum(s);
end
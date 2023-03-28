function [P]=siteprobdistribution(sequence,h,e,ms,type)
    %ms=mutation site
    [q,N]=size(h);
    P=zeros(q,1);
    
    for s=1:(ms-1)
        P=P+e((s-1)*q+sequence(s),(ms-1)*q+(1:q))';
    end
    
    P=P+h(:,ms);
    
    for s=(ms+1):N
        P=P+e((ms-1)*q+(1:q),(s-1)*q+sequence(s));
    end
    P=exp(P);
    P=P/sum(P);
    
    if (type==2)
       P=cumsum(P);
    end
    
end

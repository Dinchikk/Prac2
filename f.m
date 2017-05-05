function [ res] = f(psi , alpha )
    if((psi(2) < alpha) && (psi(2) > 0))
        res = psi(2);
    else
        if(psi(2)< 0)
            res = 0;
        else
            res = alpha;
        end
        
    end       


end


function [ res] = u1_0(t,y,alpha)
    
   
    if(y(4) < 0)
        res = 0;
    else
        if(y(4)<alpha)
            res = y(4);
        else
            res = alpha;
        end
    end
    
    %res = y(4);
end


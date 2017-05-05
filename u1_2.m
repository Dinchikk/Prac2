function [ res] = u1_2(t,y)
    
    alpha = 0.1;
    if(y(4) < 0)
        res = alpha;
    else
        res = 0;
    end
end


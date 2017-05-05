function [ res] = u2_0(t,y)
    k1 = 8;
    k2 = 9;
    if(y(2)*y(4)) < 0
        res = k2;
    else
        res = k1;
    end
end


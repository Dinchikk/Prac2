function [ res] = u2_1(t,y,k1 , k2)
    %k1 = 7;
    %k2 = 8;
    if(y(2)*y(4)) < 0
        res = k2;
    else
        res = k1;
    end
end



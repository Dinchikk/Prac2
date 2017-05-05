function res = u2(t_star , t , k1 , k2 )
    if(t < t_star)
        res = k1;
    else
        res = k2;
    end

end


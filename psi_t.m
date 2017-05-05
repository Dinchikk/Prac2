function [ res ] = psi_t(A , t , psi0 )
    
    [~ , psi] = ode45(@(tau, x) A(tau) * x , linspace(0, t, 100), psi0);
    
    x1_end = psi(end , 1);
    x2_end = psi(end , 2);
    
    if(abs(x1_end) > 5 )
        if(x1_end < 0)
            x1_end = -5;
        else
            x1_end = 5;
        end
    end
    
    if(abs(x2_end) > 5 )
        if(x2_end < 0)
            x2_end = -5;
        else
            x2_end = 5;
        end
    end
    
    res = [x1_end , x2_end];
end


function res = func(psi0 , T , k1 , p , x2_end , L)
    res = zeros(1 , 2);
    g = 9.8;
    
    res(1) = exp(-0.5*T*(1+k1))/(2*(1+k1)*p)*( (exp(T*(1+k1)) - 1)*p*psi0(2)*cosh(0.5*p*T) - (2*psi0(1)*( exp(T*(1+k1)) - 1 ) - (1+k1)*(exp(T*(1+k1)) + 1 )*psi0(2) )* sinh(0.5*p*T) ) - x2_end;    
    res(2) = exp(0.5*T*(1 + k1))/p*( (-2*psi0(1) + (1 + k1)*psi0(2))*sinh(0.5*p*T) + p*psi0(2)*cosh(0.5*p*T));
end


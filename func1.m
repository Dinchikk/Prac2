function res = func1(psi0 , t_star,t_alpha , k1 , p ,alpha)
    res = zeros(1 , 2);
    
    res(1) = ( ((1 + k1)*sinh(0.5*p*t_alpha) + p*cosh(0.5*p*t_alpha) )*psi0(2) - alpha*p*exp(-0.5*t_alpha*(1 + k1)))/(2*sinh(0.5*p*t_alpha)) - psi0(1);
    res(2) = tanh(0.5*p*t_star)  - p*psi0(2)/(2*psi0(1) - (1 + k1)*psi0(2));
end


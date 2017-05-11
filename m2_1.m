%%




clear
T = 2.4;
eps = 0.2;
L = -1.6;
g = 8;
k1 = 7;
k2 = 9;
alpha = 0.1;

J_min = 100000;



p = (k1*k1 + 2*k1 -4*g + 1);
if(p < 0)
    error('err:  p < 0')
else
   p = p^0.5;
end;

hold on


%% 

n = 500;

x2_end = linspace(0 , eps , n);


plot(linspace(-eps , eps,100) , eps*ones(1 , 100), 'g');
plot(linspace(-eps , eps,100) , -eps*ones(1 , 100), 'g');
plot(eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');
plot(-eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');
for(i = [1:n])
    
     fun = @(x)func(x , T , k1 , p , x2_end(i) , L);
     x0 = [0,0];         
     psi0 = fsolve(fun,x0);
     
     A = @(tau, x)[0 1 0 0; -g -(1+u2_1(tau,x,k1 , k2)) 0 0;  0 0 0 g; 0 0 -1 1+u2_1(tau,x,k1 , k2) ];
     f = @(tau, x) [0 ; u1_1(tau , x);0;0];

     [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi0(1);psi0(2)]);
    
     plot(X(:,1), X(:,2));
     %disp(X(end , 1));
     %disp(X(end , 2));
     
     
     
    if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
        u1_2 = zeros(1 , max(size(X)));
        u1_2= X(:,4);
        u1_2 = u1_2.^2;
       % disp(u1_2);
       J_tmp = trapz(t_time , u1_2);
       
       %disp(J_tmp);
       if (J_tmp < J_min)
           J_min = J_tmp;
           x_min = X;
           time_min = t_time;
           plot(X(:,1) , X(: , 2) , 'y');
           u1_min = X(:,4);
           u2_min = ones(1 , max(size(X)));
           for(i2 = [1:max(size(X))])
              u2_min(i2) = u2_1(t_time(i2),X(i2 , :),k1 , k2); 
           end
       end
           
    end 
      
    

end



ylabel('$$x_2$$','interpreter','latex','fontsize',15,'rotation',0);
xlabel('$$x_1$$','interpreter','latex','fontsize',15);
hold off

if(J_min < 10000)
    figure(2)
    hold on
    plot(x_min(: , 1) ,x_min(: , 2) , 'b' );
    plot(linspace(-eps , eps,100) , eps*ones(1 , 100), 'g');
    plot(linspace(-eps , eps,100) , -eps*ones(1 , 100), 'g');
    plot(eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');
    plot(-eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');
    ylabel('$$x_2$$','interpreter','latex','fontsize',13,'rotation',0);
    xlabel('$$x_1$$','interpreter','latex','fontsize',13);
    hold off


    figure(3)
    hold on
    plot(time_min ,x_min(: , 1) , 'b' );
    ylabel('$$x_1$$','interpreter','latex','fontsize',13,'rotation',0);
    xlabel('$$t$$','interpreter','latex','fontsize',13);
    hold off


    figure(4)
    hold on
    plot(time_min ,x_min(: , 2) , 'b' );
    ylabel('$$x_2$$','interpreter','latex','fontsize',13,'rotation',0);
    xlabel('$$t$$','interpreter','latex','fontsize',13);
    hold off
    disp(J_min);

    figure(5)
    hold on
    plot(time_min ,u1_min , 'b' );
    ylabel('$$u_1$$','interpreter','latex','fontsize',13,'rotation',0);
    xlabel('$$t$$','interpreter','latex','fontsize',13);
    hold off
    
    figure(6)
    hold on
    plot(time_min ,u2_min , 'b' );
    ylabel('$$u_2$$','interpreter','latex','fontsize',13,'rotation',0);
    xlabel('$$t$$','interpreter','latex','fontsize',13);
    hold off
    
    
    disp('Минимизатор:')
    disp(J_min);

else
    disp('Нет решения');
end

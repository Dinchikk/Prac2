%%

clear
T = 2;
eps = 0.2;
L = -1;
g = 9.8;
k1 = 8;
k2 = 9;
alpha = 0.1;

J_min = 100000;

n = 20;

p = (k1*k1 + 2*k1 -4*g + 1);
if(p < 0)
    disp('err');
else
   p = p^0.5;
end;

hold on


plot(linspace(-eps , eps,100) , eps*ones(1 , 100), 'g');
plot(linspace(-eps , eps,100) , -eps*ones(1 , 100), 'g');
plot(eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');
plot(-eps*ones(1 , 100),linspace(-eps , eps,100) , 'g');



%% 
n = 100;
A = [0 1; -g -(1 + k1)];
f = [0 ; alpha];
time = linspace(0 , T ,n);


[t_tmp , x] = ode45(@(tau, x) A * x + f, time, [L ; 0]);
plot(x(:,1) , x(: , 2) , 'b');



if ( ( abs(x(end , 1)) < eps) && (abs(x(end , 2)) < eps) )
    u1_2 = alpha^2 * ones(1 , n);
    J_tmp = trapz(time ,u1_2);
    if(J_tmp < J_min)
       J_min = J_tmp;
       x_min = x;
       time_min = t_tmp;
       u1_min = ones(1 , max(size(t_tmp))).*alpha;
       u2_min = ones(1 , max(size(t_tmp))).*k1;
    end
end


%% 


n = 8;
psi2_0 = linspace(0 , alpha , n);
t_alpha = linspace(0 , T , n);


for(i = [2:1:(n-1)])
   for(j = [1:1:(n-1)])
       
       psi1_0 = 1/(2*sinh(0.5*p*t_alpha(i)))*(( (1 + k1)*sinh(0.5*p*t_alpha(i)) + p.*cosh(0.5*p*t_alpha(i))).*psi2_0(j) - alpha*p*exp(-0.5*t_alpha(i)*(1 + k1)));
       
       
     A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
     f = @(tau, x) [0 ; u1_0(tau , x);0;0];

     [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0(j)]);
    
     plot(X(:,1), X(:,2),'m');

     
    if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps)) 
        u1_2 = ones(1 , max(size(X)));
        for(i2 = [1:max(size(X))])
          u1_2(i2) = u1_0(t_time(i),X(i , :)); 
        end
        J_tmp = trapz(t_time , u1_2.^2);
       %disp(J_tmp);
       if (J_tmp < J_min)
           J_min = J_tmp;
           x_min = X;
           time_min = t_time;
           plot(X(:,1) , X(: , 2) , 'y');
           u1_min =  u1_2;
           u2_min = ones(1 , max(size(X)));
           for(i2 = [1:max(size(X))])
              u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
           end
       end
           
    end 
      
   end
       
    
end

%% 

n = 8;
t_star = linspace(0.01 , T , n);
t_alpha = linspace(0.01 , T , n);

for(i = [1:(n-1)])
    for(j = [(i+1):n])
        psi2_0 = alpha * exp(-0.5*t_alpha(i)*(1 + k1))*sinh(0.5*p*t_star(j))/sinh(0.5*p*(t_star(j) - t_alpha(i)));
        psi1_0 = 1./(2*sinh(0.5*p*t_alpha(i))).*(( (1 + k1)*sinh(0.5*p*t_alpha(i)) + p.*cosh(0.5*p*t_alpha(i))).*psi2_0 - alpha*p*exp(-0.5*t_alpha(i)*(1 + k1)));
        
        A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
     f = @(tau, x) [0 ; u1_0(tau , x);0;0];

     [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0]);
    
     plot(X(:,1), X(:,2),'b');
     %disp(X(end , 1));
     %disp(X(end , 2));
     
     
     
    if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
         u1_2 = ones(1 , max(size(X)));
        for(i2 = [1:max(size(X))])
          u1_2(i2) = u1_0(t_time(i),X(i , :)); 
        end
        J_tmp = trapz(t_time , u1_2.^2);
       %disp(J_tmp);
       if (J_tmp < J_min)
           J_min = J_tmp;
           x_min = X;
           time_min = t_time;
           plot(X(:,1) , X(: , 2) , 'y');
           u1_min =  u1_2;
           u2_min = ones(1 , max(size(X)));
           for(i2 = [1:max(size(X))])
              u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
           end
       end
           
    end 
      

    end
end


%%

n = 8;
t_star = linspace(0.0001 , T , n);
psi2_0 = linspace(0.0001 , alpha , n);


for(i = [1:n])
    for(j = [1 : n])
        psi1_0 = 0.5*psi2_0(i)*(p*cosh(0.5*p*t_star(j)) + 1 +k1);
        
        A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
        f = @(tau, x) [0 ; u1_0(tau , x);0;0];

        [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0(j)]);
    
        plot(X(:,1), X(:,2));
        %disp(X(end , 1));
        %disp(X(end , 2));
     
     
     
        if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
             u1_2 = ones(1 , max(size(X)));
            for(i2 = [1:max(size(X))])
                u1_2(i2) = u1_0(t_time(i),X(i , :)); 
            end
            J_tmp = trapz(t_time , u1_2.^2);
       
           %disp(J_tmp);
           if (J_tmp < J_min)
               J_min = J_tmp;
               x_min = X;
               time_min = t_time;
               plot(X(:,1) , X(: , 2) , 'y');
               u1_min =  u1_2;
               u2_min = ones(1 , max(size(X)));
               for(i2 = [1:max(size(X))])
                  u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
               end
           end

        end 

    end
end    



%% 


n = 8;
t_alpha = linspace(0.0001 , T , n);
psi2_0 = linspace(0.001 , alpha , n);

for(i = [1:n])
    for(j = [1 : n])
        psi1_0 = 1/(2*sinh(0.5*p*t_alpha(j)))*(((1+k1)*sinh(0.5*p*t_alpha(j))+ p*cosh(0.5*p*t_alpha(j)) )*psi2_0(i) - alpha*p*exp(-0.5*t_alpha(j)*(1+k1)));
        
        A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
        f = @(tau, x) [0 ; u1_0(tau , x);0;0];

        [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0(i)]);
    
        plot(X(:,1), X(:,2));
        %disp(X(end , 1));
        %disp(X(end , 2));
     
     
     
        if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
             u1_2 = ones(1 , max(size(X)));
            for(i2 = [1:max(size(X))])
                u1_2(i2) = u1_0(t_time(i),X(i , :)); 
            end
            J_tmp = trapz(t_time , u1_2.^2);
       
           %disp(J_tmp);
           if (J_tmp < J_min)
               J_min = J_tmp;
               x_min = X;
               time_min = t_time;
               plot(X(:,1) , X(: , 2) , 'y');
               u1_min =  u1_2;
               u2_min = ones(1 , max(size(X)));
               for(i2 = [1:max(size(X))])
                  u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
               end
           end

        end
    end
end    



%%

n = 100;
A = [0 1; -g -(1 + k1)];
f = [0 ; 0];
time = linspace(0 , T ,n);


[t_tmp , x] = ode45(@(tau, x) A * x + f, time, [L ; 0]);

if ( ( abs(x(end , 1)) < eps) && (abs(x(end , 2)) < eps) )
    u1_min = ones(1 , max(size(t_tmp))).*0;
    u2_min = ones(1 , max(size(t_tmp))).*k1;
    u1_2 = 0 * ones(1 , n);
    J_tmp = trapz(time ,u1_2.^2);
    if(J_tmp < J_min)
       J_min = J_tmp;
       x_min = x;
       time_min = t_tmp;
       plot(x(:,1) , x(: , 2) , 'b');
    end
end


%% 

n = 8;
t_star = linspace(0.01 , T , n);
psi2_0 = 1;

for(i = [1:n])
        
        psi1_0 = 0.5*psi2_0*(p*cosh(0.5*p*t_star(i)) + 1 + k1);
        u1_psi0(7 , [1 1 1 1]);
        
        A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
        f = @(tau, x) [0 ; u1_psi0(tau , x);0;0];

        [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0]);
    
        plot(X(:,1), X(:,2),'b');
        %disp(X(end , 1));
        %disp(X(end , 2));
     
     
     
        if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
             u1_2 = ones(1 , max(size(X)));
            for(i2 = [1:max(size(X))])
                u1_2(i2) = u1_0(t_time(i),X(i , :)); 
            end
            J_tmp = trapz(t_time , u1_2.^2);
       
           %disp(J_tmp);
           if (J_tmp < J_min)
               J_min = J_tmp;
               x_min = X;
               time_min = t_time;
               plot(X(:,1) , X(: , 2) , 'y');
               u1_min =  u1_2;
               u2_min = ones(1 , max(size(X)));
               for(i2 = [1:max(size(X))])
                  u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
               end
           end

        end
end  


%% 

n = 8;
t_star = linspace(0.01 , T , n);
psi2_0 = -1;

for(i = [1:n])
        
        psi1_0 = 0.5*psi2_0*(p*cosh(0.5*p*t_star(i)) + 1 + k1);
        
        
        A = @(tau, x)[0 1 0 0; -g -(1+u2_0(tau,x)) 0 0;  0 0 0 g; 0 0 -1 1+u2_0(tau,x) ];
        f = @(tau, x) [0 ; u1_psi0(tau , x);0;0];

        [t_time , X] = ode45(@(tau, x) A(tau , x) * x + f(tau , x) , linspace(0, T, 100), [L ; 0;psi1_0;psi2_0]);
    
        plot(X(:,1), X(:,2), 'b');
        %disp(X(end , 1));
        %disp(X(end , 2));
     
     
     
        if((abs(X(end , 1)) < eps) && (abs(X(end , 2)) < eps))
             u1_2 = ones(1 , max(size(X)));
            for(i2 = [1:max(size(X))])
                u1_2(i2) = u1_0(t_time(i),X(i , :)); 
            end
            J_tmp = trapz(t_time , u1_2.^2);
       
           %disp(J_tmp);
           if (J_tmp < J_min)
               J_min = J_tmp;
               x_min = X;
               time_min = t_time;
               plot(X(:,1) , X(: , 2) , 'y');
               u1_min =  u1_2;
               u2_min = ones(1 , max(size(X)));
               for(i2 = [1:max(size(X))])
                  u2_min(i2) = u2_1(t_time(i),X(i , :),k1 , k2); 
               end
           end

        end
end  

%% 

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

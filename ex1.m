%% Tests
close all;
% integral between 0 and PI/2 of sin(x)dx = 1
% sin_x = @(x) sin(x);

sin_x = @(x) sin(x);
x_0 = 0;
x_final = pi/2;
nb_iteration = 4;

nb_iteration_convergence = 8;


% Quad midpoint
integral = quad_midpoint(sin_x, x_0, x_final, nb_iteration);
disp("Calcul with midpoint rule and 4 iterations");
disp(integral);

% Quad trapezoidal
integral = quad_trapezoidal(sin_x, x_0, x_final, nb_iteration);
disp("Calcul with trapezoidal rule and 4 iterations");
disp(integral);

% Quad simpsons
integral = quad_simpsons(sin_x, x_0, x_final, nb_iteration);
disp("Calcul with simpsons rule and 4 iterations");
disp(integral);


% Convergence rate

method1 = @(f,a,b,N) quad_midpoint(f,a,b,N);
method2 = @(f,a,b,N) quad_trapezoidal(f,a,b,N);
method3 = @(f,a,b,N) quad_simpsons(f,a,b,N);

figure;
subplot(3,1,1);
plot_convergence_error(sin_x, x_0,x_final,nb_iteration_convergence,method1);
title('Convergence error: midpoint')
xlabel('n') 
ylabel('midpoint convergence rate') 

subplot(3,1,2);
plot_convergence_error(sin_x, x_0,x_final,nb_iteration_convergence,method2);
title('Convergence error: trapezoidal')
xlabel('n') 
ylabel('trapezoidal convergence rate') 


subplot(3,1,3);
plot_convergence_error(sin_x, x_0,x_final,nb_iteration_convergence,method3);
title('Convergence error: simpsons')
xlabel('n') 
ylabel('simpsons convergence rate')


% Resolution of differential equations - tests

% Q1 
figure;
subplot(3,1,1);

y0 = 0.1;
T = 4;
N = 25;
real_y = @(t, y0) (y0 ./ sqrt(y0^2 - (y0^2 - 1) * exp(-2*t)));
f = @(x) (x-x^3);

time = linspace(0,T,N);
results_calculated = euler_method(f,y0, T, N);
results_real = real_y(time,y0);

plot(time,results_calculated,'+',time,results_real,'r-');
title('Forward Euler method with N=25, T=4 and y0=0.1')
xlabel('t') 
ylabel('y(t)')
legend({'Estimated','Real'},'Location','southeast')
legend('boxoff')

subplot(3,1,2);
results_calculated = backward_euler_method(f,y0, T, N);
plot(time,results_calculated,'+',time,results_real,'r-');
title('Backward Euler method with N=25, T=4 and y0=0.1')
xlabel('t') 
ylabel('y(t)')
legend({'Estimated','Real'},'Location','southeast')
legend('boxoff')

subplot(3,1,3);
results_calculated = crank_nicolson_method(f,y0, T, N);
plot(time,results_calculated,'+',time,results_real,'r-');
title('Crank-Nicolson method with N=25, T=4 and y0=0.1')
xlabel('t') 
ylabel('y(t)')
legend({'Estimated','Real'},'Location','southeast')
legend('boxoff')


% Q2
T=1;
p=0.01;
y0=0.1;
N = 20;

order_euler  = calculate_order(f,y0,T,N,p, real_y, @euler_method);
disp('Order of forward euler');
disp(order_euler);

order_backward_euler = calculate_order(f,y0,T,N,p, real_y, @backward_euler_method);
disp('Order of backward euler');
disp(order_backward_euler);

order_crank_nicolson = calculate_order(f,y0,T,N,p, real_y, @crank_nicolson_method);
disp('Order of Crank-Nicolson method');
disp(order_crank_nicolson);



%% Functions
function [sum] = quad_midpoint(f, a, b, N)
    h = (b - a) / N;
    % Calculation of extremity points
    sum = f(a)*h/2 + f((N+1/2)*h)*h/2;
    
    % Calculation of intermediate points
    for i = 2:N
        x_i = (i-(1/2))*h;
        sum = sum + (f(x_i)*h);
    end
end

function [sum] = quad_trapezoidal(f, a, b, N)
    h = (b - a) / N;
    % Calculation of extremity points
    sum = f(a)*h/2 + f(N*h)*h/2;
    
    % Calculation of intermediate points
    for i = 2:N
        x_i = (i-1)*h;
        sum = sum + (f(x_i)*h);
    end
end


function [sum] = quad_simpsons(f, a, b, N)
    h = (b - a) / N;
    % Calculation of extremity points
    sum = f(a)*h/6 + f(N*h)*h/6;
    
    % Calculation of intermediate points
    for i = 2:(2*N)
        x_i = (i-1)*h/2;
        if mod(i,2) == 0
            sum = sum + (f(x_i)*(4*h/6));
        else
            sum = sum + (f(x_i)*(2*h/6));
        end
    end
end

function [convergence_rate] = calculate_convergence_rate_q1(f, a, b, N, method)
    real_integral_h = integral(f,a,b);

    calculated_integral_h = method(f,a,b,N);
    error_h = abs(real_integral_h - calculated_integral_h);
    
    calculated_integral_h_2 = method(f,a,b,2*N);
    error_h_2 = abs(real_integral_h - calculated_integral_h_2);

    convergence_rate = log2(error_h/error_h_2);
end

function plot_convergence_error(f, a, b, N, method)
    n = 1;
    convergence_rates = zeros(1,N);
    n_values = zeros(1,N);
    for k = 1:N
        convergence_rates(k) = calculate_convergence_rate_q1(f, a, b, n, method);
        n_values(k) = n;
        n = n*2;
    end
    loglog(n_values,convergence_rates)
    grid on
end

%% Resolution of differential equations - functions
% y0 => result of function at 0
% N => number of steps
% T => final time
function [y] = euler_method(f,y0,T,N)
    % Method implemented : y_n_1 = yn + h*f(tn,yn)
    y = zeros(1,N);
    y(1) = y0;
    % delta time
    h = T/N;
    for k = 2:N
        y(k) = y(k-1) + h * f(y(k-1));
    end
end

function [y] = backward_euler_method(f,y0,T,N)
    % Method implemented : y_n_1 = yn + h*f(tn_1,yn_1)
    y = zeros(1,N);
    h = T/N;
    y(1) = y0;
    for i = 1:N-1
        y(i+1) = secant(@(Y) y(i) + h*f(Y) - Y, i*h, i*(h+1), 1e-06);
    end
end

function [y] = crank_nicolson_method(f,y0,T,N)
    % Method implemented : y_n_1 = yn + 0.5*h*f(tn,yn) + 0.5*h*f(tn_1,yn_1)
    y = zeros(1,N);
    h = T/N;
    y(1) = y0;
    for i = 1:N-1
        y(i+1) = secant(@(Y) y(i) + 0.5*h*f(y(i)) + 0.5*h*f(Y) - Y, i*h, i*(h+1), 1e-06);
    end
end


function [x_k_p_1] = secant(f,x0,x1,precision)
    zero_condition = 1e-10;
    x_k_m_1 = x0;
    x_k = x1;
    x_k_p_1 = x_k - f(x_k) * (x_k-x_k_m_1) / (f(x_k)-f(x_k_m_1));
    while abs(f(x_k_p_1)) > precision
        x_k_m_1 = x_k;
        x_k = x_k_p_1;
        if abs(f(x_k)-f(x_k_m_1)) <= zero_condition
            error('division by zero, cannot find a solution');
        else
            x_k_p_1 = x_k - f(x_k) * (x_k-x_k_m_1) / (f(x_k)-f(x_k_m_1));
        end
    end
end

%% Asymptotic order
function [changement,order] = check_changement(previous_rate, new_rate, precision)
    changement = true;
    if (abs(previous_rate - new_rate) <= precision) || (abs(new_rate - round(new_rate)) <= precision)
        changement = false;
    end 
    order = round(new_rate);
end

function [convergence_rate] = calculate_convergence_rate_q2(method,y0,T,N,real_y,f)
    calculated_Y = method(f,y0,T,N);
    calculated_y_2N = method(f,y0,T,2*N);
    error_1 = abs(real_y - calculated_Y(end));
    error_2 = abs(real_y - calculated_y_2N(end));
    convergence_rate = log2(error_1/error_2);
end


function [order] = calculate_order(f,y0,T,N,p, real_y, method)
    changement = true;
    convergence_rates = zeros(1,N);
    for i=1:N
        convergence_rates(i) = calculate_convergence_rate_q2(method,y0,T,i,real_y(T,y0),f);
        if(i > 1)
            [changement,order] = check_changement(convergence_rates(i-1),convergence_rates(i),p);
        end
        if(~changement)
            break;
        end
    end
end


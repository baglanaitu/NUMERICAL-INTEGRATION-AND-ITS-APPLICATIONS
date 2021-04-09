% Intialization of the variables
T = 1; % time value of 1 as given in the question
y0 = 0.1; % intial value for y
N = 1:100; % Possible N values to make h approaches 0 
p = 0.01; % percision value
r_eular = []; % keep track of convergence rate of Methods
r_backward = [];
r_crank = [];
y = @(t) y0./(sqrt(y0^2-(y0^2-1).*exp(-2*t))); % the anatytical function
h_eular = 0; % keep track of the current N for h goes to 0
h_backward = 0;
h_crank = 0;
o_eular = 0; % Keep track of Method order
o_backward = 0;
o_crank = 0;
c_eular = true; % Keep track of method changing status
c_backward = true;
c_crank = true;
% looking for convergence rate of Eular Method
for i = N
    r_eular = [r_eular  eval_conver_rate(@eular_method,y0,T,i,y(1))];
    if(i > 1)
        [c_eular,o_eular] = isChanging(r_eular(i-1),r_eular(i),p);
    end
    if(~c_eular)
        h_eular = i;
        break;
    end
end
% looking for convergence rate of Backward Eular Method
for i = N
    r_backward = [r_backward  eval_conver_rate(@backward_eular_method,y0,T,i,y(1))];
    if(i > 1)
        [c_backward,o_backward] = isChanging(r_backward(i-1),r_backward(i),p);
    end
    if(~c_backward)
        h_backward = i;
        break;
    end
end
% looking for convergence rate of Crank-Nicolson Method
for i = N
    r_crank = [r_crank  eval_conver_rate(@crank_nicolson_method,y0,T,i,y(1))];
    if(i > 1)
        [c_crank,o_crank] = isChanging(r_crank(i-1),r_crank(i),p);
    end
    if(~c_crank)
        h_crank = i;
        break;
    end
end

%Displaying the Assymptotic Order along with number of h
disp("The assumptotic order of Eular Method:");
disp(o_eular);
disp("The assumptotic order of Backward Eular Method:");
disp(o_backward);
disp("The assumptotic order of Crank Nicolson Method:");
disp(o_crank);

% Plateau Checker
function [c,o] = isChanging(r_prev,r_new,p)
    c = true;
    o = inf;
    if (abs(r_prev - r_new) <= p)
        c = false;
    elseif(abs(r_new - round(r_new)) <= p)
        c = false;
    end 
    o = round(r_new);
end
%Eular Method
function [y,t] = eular_method(y0,T,N)
h = T / N;
yn = zeros(N+1,1);
tn = zeros(N+1,1);
yn(1) = y0;
% assuming time start by 0
tn(1) = 0;
for i = 1:length(yn)-1
    yn(i+1) = yn(i) + h *(yn(i) - yn(i)^3);
    tn(i+1) = tn(i) + h;
end
 y = yn;
 t = tn;
end

%backward Eular Method
% To solve the nonlinear eq., we will use the secant method
function [y,t] = backward_eular_method(y0,T,N)
h = T / N;
yn = zeros(N+1,1);
tn = zeros(N+1,1);
yn(1) = y0;
% assuming time start by 0
tn(1) = 0;
for i = 1:length(yn)-1
    tn(i+1) = tn(i) + h;
    [yn(i+1),~] = secant(@(x) yn(i) + h *(x - x^3) - x,tn(i),tn(i+1), 0.00001);
end
 y = yn;
 t = tn;
end


% Crank-Nicolson Method
function [y,t] = crank_nicolson_method(y0,T,N)
h = T / N;
yn = zeros(N+1,1);
tn = zeros(N+1,1);
yn(1) = y0;
% assuming time start by 0
tn(1) = 0;
for i = 1:length(yn)-1
    tn(i+1) = tn(i) + h;
    [yn(i+1),~] = secant(@(x) yn(i) + 0.5 * h * (yn(i) - yn(i)^3)  + 0.5 * h * (x - x^3) - x,tn(i),tn(i+1), 0.00001);
end
 y = yn;
 t = tn;
end

%Emplementing asymptotic order evaluator
function [r] = eval_conver_rate(f,y0,T,N,y_exact)
    [y_aprrox,~] = f(y0,T,N);
    [y_aprrox_2,~] = f(y0,T,2*N);
     err_1 = abs(y_exact - y_aprrox(end));
     err_2 = abs(y_exact - y_aprrox_2(end));
     r = log2(err_1/err_2);
end




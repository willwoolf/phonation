%% Two Mass Model - Normalised ODE Attempt 3

%final nondimensionalisation, HOPEFULLY

%set the required constants

% alpha is quotient of second mass by first mass, we specify only its
% inverse alpha_inv = m1/m2 since we only ever divide by alpha.
% lambda is the quotient of their support spring stiffnesses k2/k1.
% beta is proportional to the forcing pressure at lung.
% omega is the quotient of the coupling stiffness by spring 1, i.e. if we
% know omega and lambda, you can find the stiffness k_2.

alpha_inv = 1;
lambda = 0.8;
beta = 3;
omega = 0.3;

%% Stationary points and general equilibria solvers

%one-dimensional solvers for individual stationary points, neglecting the
%stiffness coupling
stat_point_11 = fzero(@(x) 1 - x + beta*(1 - power(x, -2)),[1e-5 power(2*beta, 1/3)]);
stat_point_12 = fzero(@(x) 1 - x + beta*(1 - power(x, -2)),[power(2*beta, 1/3), 10]);
stat_point_21 = fzero(@(x) lambda*(1 - x) + beta*(1 - power(x, -2)),[1e-5 power(2*beta/lambda, 1/3)]); 
stat_point_22 = fzero(@(x) lambda*(1 - x) + beta*(1 - power(x, -2)),[power(2*beta, 1/3), 10]);

%solver problem to find equilibria. MUST take vertical vector input
init_search_point = [stat_point_12; stat_point_22];
[steady_soln,fval] = fsolve(@(v) zeroesFunction2(v,lambda,beta,omega),init_search_point);

x_equilibrium = steady_soln(1);
y_equilibrium = steady_soln(2);

%F = zeroesFunction(x_equilibrium,y_equilibrium,lambda,beta,omega)
%J = zeroesGradient(x_equilibrium,y_equilibrium,lambda,beta,omega)

%% Bifurcation Diagram

beta_arr = linspace(0.3,0.6,21)
lambda;
omega;

%use previous computations
solns_x = double.empty(2,2,0)
solns_y = double.empty(2,2,0)
for beta = beta_arr
    starting_point_mat = findStatPoints(lambda,beta);
    stat_point_11 = starting_point_mat(1,1);
    stat_point_12 = starting_point_mat(1,2);
    stat_point_21 = starting_point_mat(2,1);
    stat_point_22 = starting_point_mat(2,2);
    [kilib_x, kilib_y] = findEquilibria(stat_point_11,stat_point_12,stat_point_21,stat_point_22,lambda,beta,omega);
    solns_x = cat(3,solns_x,kilib_x);
    solns_y = cat(3,solns_y,kilib_y);
end
equilibrium_11_x = reshape(solns_x(1,1,:),1,[]);
equilibrium_12_x = reshape(solns_x(1,2,:),1,[]);
equilibrium_21_x = reshape(solns_x(2,1,:),1,[]);
equilibrium_22_x = reshape(solns_x(2,2,:),1,[]);
equilibrium_11_y = reshape(solns_y(1,1,:),1,[]);
equilibrium_12_y = reshape(solns_y(1,2,:),1,[]);
equilibrium_21_y = reshape(solns_y(2,1,:),1,[]);
equilibrium_22_y = reshape(solns_y(2,2,:),1,[]);

hold on
plot(beta_arr,equilibrium_11_x)
plot(beta_arr,equilibrium_12_x)
plot(beta_arr,equilibrium_21_x)
plot(beta_arr,equilibrium_22_x)
hold off


%% Eigenvalue computation

df1dx = partDf1(x_equilibrium,beta,omega);
df1dy = omega;
df2dx = omega;
df2dy = partDf2(y_equilibrium,lambda,beta,omega);

jcbn = [
    0, 1, 0, 0;
    df1dx, 0, df1dy, 0;
    0, 0, 0, 1;
    alpha_inv*df2dx, 0, alpha_inv*df2dy, 0;
];

eig_vals = eig(jcbn)

hold on
for val = eig_vals
    plot([-5 5],[0 0],'k-')
    plot([0 0],[-5 5],'k-')
    plot(real(val),imag(val),'rx','MarkerSize',8)
end
hold off

%% Eigenvalue graphics (computation +)

% plots the eigenvalues of a general equilibrium solution as a parameter
% changes. the variable of interest must be substituted into the matrix.

% the implementation of newtonSolver needs to be moved, the below code is
% finding the jacobian about a fixed point in space for different
% parameters, not accounting for changing locations of the equilibrium
% solutions as the parameters change in the first place

init_search_point = [stat_point_12;stat_point_22];
eigenvals = double.empty(4,0);
interesting_params = [0.25 0.5 1 2 4];
for param = interesting_params
    [steady_soln,fval] = fsolve(@(v) zeroesFunction2(v,param,beta,omega),init_search_point);
    x_steady = steady_soln(1)
    y_steady = steady_soln(2)
    jcb = [
        0, 1, 0, 0;
        partDf1(x_steady,beta,omega), 0, omega, 0;
        0, 0, 0, 1;
        omega, 0, partDf2(y_steady,param,beta,omega), 0
    ];
    jcb_eigs = eig(jcb);
    jcb_eigs = sort(jcb_eigs)
    eigenvals = [eigenvals jcb_eigs];
end

eigenreal = real(eigenvals);
eigenimag = imag(eigenvals);
hold on
subplot(121)
plot(eigenreal(1,:),eigenimag(1,:),'rs--')
plot(eigenreal(2,:),eigenimag(2,:),'bo--')
subplot(122)
plot(eigenreal(3,:),eigenimag(3,:),'rs--')
plot(eigenreal(4,:),eigenimag(4,:),'bo--')
hold off

%% Euclidean (2-norm) distance of the computation from the unstable equilibrium (1,1).

%4 dimensional norm from the equilibrium
unsteady = [1,0,1,0]';
X=x';
r = vecnorm(X-unsteady);

%2 dimensional norm from the position
position = [1,1]'
Y=[x(:,1),x(:,3)]'
d = vecnorm(Y-position)


%% Contour graphic of equilibrium solutions

x_axis = [-5:0.05:15];
y_axis = [-5:0.05:15]';
z_axis_1 = 1 - x_axis + beta*(1-1./(x_axis.^2)) + omega*(y_axis-x_axis);
z_axis_2 = lambda*(1 - y_axis) + beta*(1-1./(y_axis.^2)) + omega*(x_axis-y_axis);

hold on
contour(x_axis,y_axis,z_axis_1,[-0.001,0,0.001],'b-')
contour(x_axis,y_axis,z_axis_2,[-0.001,0,0.001],'r-')
plot([-5 15],[-5 15],'g-')
hold off

%% Implementation of ODE45

%enter initial conditions
opts = odeset('AbsTol',1e-20);
tspan = [0 100];
init = [stat_point_12, 0, stat_point_22, 0];
[t,x] = ode45(@(t,x) twoMassModel(t,x,alpha_inv, lambda, beta, omega),tspan,init,opts);

E = totalEnergy(x(:,1),x(:,2),x(:,3),x(:,4),alpha_inv,lambda,beta,omega);

hold on
%plot(x(:,1), x(:,2));
%plot(x(:,3), x(:,4));

plot(t, x(:,1))
%plot(t, x(:,3))

%plot(t,E)

%plot(x(:,1),x(:,3))
hold off

%% Les functiones

%the ODE function on two masses taking parameters alpha, lambda, beta,
%omega
function dhdt = twoMassModel(t,x,alpha_inv,lambda,beta,omega)

%take inputs, initialise the h variables
dhdt = zeros(size(x));

h_1 = x(1);
dh_1 = x(2);
h_2 = x(3);
dh_2 = x(4);

%compute the equation

dhdt(1) = dh_1;
dhdt(2) = 1 - h_1 + beta*(1-1/(h_1)^2) + omega*(h_2-h_1);
dhdt(3) = dh_2;
dhdt(4) = alpha_inv*(lambda*(1 - h_2) + beta*(1-1/(h_2)^2) + omega*(h_1-h_2));

end

%vector valued function, we solve F=0 to find the equilibria
function out_vector = zeroesFunction(x,y,lambda,beta,omega)
out_x = 1-x + beta*(1-(1./(x.^2))) + omega.*(y-x);
out_y = lambda*(1-y) + beta*(1-(1./(y.^2))) + omega.*(x-y);
out_vector = [out_x; out_y];
end

function out_vector = zeroesFunction2(v,lambda,beta,omega)
x = v(1);
y = v(2);
out_x = 1-x + beta*(1-(1./(x.^2))) + omega.*(y-x);
out_y = lambda*(1-y) + beta*(1-(1./(y.^2))) + omega.*(x-y);
out_vector = [out_x; out_y];
end

%jacobian of the zeroesFunction used for Newton's method
function jacob = zeroesGradient(x,y,lambda,beta,omega)
dfdx = 2*beta*power(x,-3)- 1 - omega; %derivative of f1
dfdy = 2*beta*power(y,-3) - lambda - omega; %drv of f2
jacob = [dfdx omega; omega dfdy];
end

%partial derivatives of f1 and f2
function out = partDf1(x,beta,omega)
    out = -1 + 2*beta*power(x,-3) - omega;
end
function out = partDf2(x,lambda,beta,omega)
    out = -lambda + 2*beta*power(x,-3) - omega;
end

% conserved energy term

function vec = totalEnergy(x1,v1,x2,v2,alpha_inv,lambda,beta,omega)
kinetic_energy_term = (v1.^2 + (alpha_inv^-1)*(v2.^2))/2;
x_potential_term = (1+beta)*x1 - (x1.^2)/2 + beta./x1;
y_potential_term = (lambda+beta)*x2 - lambda*(x2.^2)/2 + beta./x2;
mixed_term = -omega*((x1-x2).^2)/2;
vec = kinetic_energy_term - x_potential_term - y_potential_term - mixed_term;
end


function starting_pts = findStatPoints(lambda,beta)
stat_point_11 = fzero(@(x) 1 - x + beta*(1 - power(x, -2)),[1e-5 power(2*beta, 1/3)]);
stat_point_12 = fzero(@(x) 1 - x + beta*(1 - power(x, -2)),[power(2*beta, 1/3), 10]);
stat_point_21 = fzero(@(x) lambda*(1 - x) + beta*(1 - power(x, -2)),[1e-5 power(2*beta/lambda, 1/3)]); 
stat_point_22 = fzero(@(x) lambda*(1 - x) + beta*(1 - power(x, -2)),[power(2*beta, 1/3), 10]);
starting_pts = [
    stat_point_11, stat_point_12;
    stat_point_21, stat_point_22
];
end

function [equilibria_x,equilibria_y] = findEquilibria(sp11,sp12,sp21,sp22,lambda,beta,omega)
%opts_2 = optimset()
steady_soln_11 = fsolve(@(v) zeroesFunction2(v,lambda,beta,omega),[sp11,sp21]);
steady_soln_12 = fsolve(@(v) zeroesFunction2(v,lambda,beta,omega),[sp11,sp22]);
steady_soln_21 = fsolve(@(v) zeroesFunction2(v,lambda,beta,omega),[sp12,sp21]);
steady_soln_22 = fsolve(@(v) zeroesFunction2(v,lambda,beta,omega),[sp12,sp22]);
if any([steady_soln_11 steady_soln_12 steady_soln_21 steady_soln_22]<0)
    equilibria_x = [
        missing missing;
        missing missing
    ];
    equilibria_y = [
        missing missing;
        missing missing
    ];
end
equilibria_x = [
    steady_soln_11(1), steady_soln_12(1);
    steady_soln_21(1), steady_soln_22(1)
];
equilibria_y = [
    steady_soln_11(2), steady_soln_12(2);
    steady_soln_21(2), steady_soln_22(2)
];
end




% sorter for eigenvalues
%function sorted = eigenSort(list_arg, memory_test)
%    arg_1 = list_arg(1)
%    arg_2 = list_arg(2)
%    arg_3 = list_arg(3)
%    arg_4 = list_arg(4)
%    %rely on size of eigenvalue
%    if abs(arg_1 - memory_test) < abs(arg_1)

%    if abs(arg_1) == abs(arg_2) && 
%        sorted = list_arg
%    elseif abs(arg_1) == abs(arg_3)
%        sorted = [arg_1 arg_3 arg_2 arg_4]
    

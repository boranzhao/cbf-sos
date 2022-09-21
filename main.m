%% synthesis of CBF using SOS and linear-like form
clear;
yalmip('clear')
% dynamics
% an unstable linear system 
% from  Clark, A., 2021, Verification and synthesis of control barrier
% functions, CDC
A = [1 0; -1 4]; B=[1;0];
% % Wang, L., Han, D. and Egerstedt, M., 2018, Permissive barrier
% % certificates for safe stabilization using sum-of-squares. ACC
% A = [0  1;-1 0]; B = [0 1];
% x = sdpvar(2,1);
n = size(A,1); nu = size(B,2);
x = sdpvar(n,1); 
x_store = x; % for later use
dynamics = @(t,x,u) A*x+B*u;

include_input_limits = 1;

% given a fixed point, rewrite the dynamics using the deviation states in linear-like form:
% fixed point
xstar = zeros(n,1);
ustar = zeros(nu,1);
% x_til = A(x_til)x_til + B(x_til)u_til, where
% x_til = x - xstar; util = u-ustar;

% x_store = x;
M = eye(n); 

% settings
deg_X = 4;
deg_Y = 2;
X_states_index = 2;

% state constraints |cx*Zx}|<1; a set C:={x||cx*Zx}|<1};
cx = x';
% input constraints |Au*u|<1;
Au = 1;

% specif a set D={x|h0(x)>=0} that contains the original safety set
h0_x = 1-(x'*x);

% polynomials and matrix coefficients
polyX = monolist(x(X_states_index),deg_X); % note that X can only depend on the states whose derivatives are not influenced by control inputs
polyY = monolist(x,deg_Y);    
d_polyX = jacobian(polyX,x(X_states_index));
X_coef = sdpvar(n,n,length(polyX),'symmetric');
Y_coef = sdpvar(nu,n,length(polyY));
X = zeros(n); Y = zeros(nu,n);
% dX_dx = zeros(n,n,length(X_states_index));
for j=1:length(X_states_index)
    dX_dx{j} = sdpvar(n,n);
end
% dX_dx = X_coef(:,:,1)*d_polyX(1,:);
% dX_dx(dX_dx-dX_dx==0)=0;
for i = 1:length(polyX)
    X = X + X_coef(:,:,i)*polyX(i);
    for j = 1:length(X_states_index)
        dX_dx{j} =  dX_dx{j} - dX_dx{j}*(i==1) + X_coef(:,:,i)*d_polyX(i,j);
    end    
end

for i = 1:length(polyY)
    Y = Y + Y_coef(:,:,i)*polyY(i);
end
paras = [X_coef(:);Y_coef(:)];

%% formulate the constraitns
Zx = x;
v = sdpvar(n,1); w = sdpvar(1);
v1 = [v;w]; 
constrs = [];

% The following constraint is redundant after having constraint (5)
% (1): positive definiteness of X in D
% epsilon = 1e-5;
% [L1,c_L1,v_L1] = polynomial([x;v],2);
% exp1 = v'*(X-epsilon*eye(n))*v - L1*h0_x;
% constrs = [sos(exp1); sos(L1)];
% paras = [paras; c_L1];

% (2): \dot h >=0 in D
X_dot = zeros(n);
for j=1:length(X_states_index)
    X_dot = X_dot + dX_dx{j}*(A(X_states_index(j),:)*Zx);    
end
tmp = M*(A*X+B*Y); 
[L2,c_L2,v_L2] = polynomial([x;v],2);
exp2 = v'*(X_dot-(tmp+tmp'))*v - L2*h0_x;
constrs = [constrs sos(L2) sos(exp2)];
paras = [paras; c_L2];

% (3): |Au*u_til|<=1;
if include_input_limits 
    exp3 = {}; L3 ={};
    for j = 1:size(Au,1)    
        Psi_j = [1  Au(j,:)*Y; (Au(j,:)*Y)' X];
        [L3{j},c_L3,v_L3] = polynomial([x;v1],2);
        exp3{j} = v1'*Psi_j*v1 - L3{j}*h0_x;   
        constrs = [constrs sos(exp3{j}) sos(L3{j})];
        paras = [paras; c_L3];
    end
end

% (4): Ch = {x|h(x)>=0} is in C:={x||cx*Zx}|<1};
exp4 = {}; L4 ={};
for i =1:size(cx,1)
    Phi_i = [1 cx(i,:)*X; (cx(i,:)*X)' X];
    [L4{i},c_L4,v_L4] = polynomial([x;v1],2);
    exp4{i} = v1'*Phi_i*v1 - L4{i}*h0_x;  
    constrs = [constrs sos(exp4{i}) sos(L4{i})];
    paras = [paras; c_L4];
end

% (5): X >= X0; 
X0 = sdpvar(n);
[L5,c_L5,~] = polynomial([x;v],2);
% obj = -logdet(X0);                        % for SDPT3: no need to inlcude X0>=0 since this automatically assumed when using logdet function, https://yalmip.github.io/command/logdet/     
obj = -geomean(X0);                         % for other solvers
exp5 = v'*(X-X0)*v-L5*h0_x;
% constrs = [constrs  sos(L5) X0<=1e-3*eye(n)];             % X0>=1e-6*eye(n)
constrs = [constrs sos(exp5) sos(L5)];      % X0>=1e-6*eye(n)
paras = [paras; c_L5];

ops = sdpsettings('solver','mosek','sos.numblkdg',1e-9); %, ,'sos.numblkdg',1e-7 sometimes this will cause infinite loop
[sol,v,Q,res] = solvesos(constrs,obj,ops,paras);
max_residual = max(res)
disp(sol.info);
%% 
if sol.problem ==0 || sol.problem == 4 
   % ------------------------- extract X_fcn  ------------------------
    X_coef = value(X_coef) 
    X_coef(abs(X_coef)<=1e-10) = 0;
    
    Y_coef = value(Y_coef);
    X_fcn = zeros(n); Y_fcn = zeros(nu,n);
    for i=1:length(polyX)
        X_fcn = X_fcn+ X_coef(:,:,i)*polyX(i);
    end
    for i=1:length(polyY)
        Y_fcn = Y_fcn+ Y_coef(:,:,i)*polyY(i);
    end
    X_fcn = clean(X_fcn, 1e-10);
    Y_fcn = clean(Y_fcn, 1e-10);
    XX = x_store; % so that the result from "sdisplay" command uses "XX"
    s = sdisplay(X_fcn);
    s2 = sdisplay(Y_fcn);
    syms XX [n 1]
    syms X_fcn [n n]
    syms Y_fcn [nu n]
    for i=1:n
        for j=1:n
            X_fcn(i,j) = eval(s{i,j});
        end
    end
    matlabFunction(X_fcn,'File','X_fcn1','Vars',{XX});
    X_fcn = matlabFunction(X_fcn,'Vars',{XX});
    
    for i=1:nu
        for j=1:n
            Y_fcn(i,j) = eval(s2{i,j});
        end
    end
    matlabFunction(Y_fcn,'File','Y_fcn1','Vars',{XX});
    Y_fcn = matlabFunction(Y_fcn,'Vars',{XX});
end
h_fcn = @(x) 1- x'/X_fcn(x)*x;
Zx_fcn = @(x) x;
return;

%% refine the CBF function
% maximize c in h = c-Zx'*X^-1*Zx such that {x|h>=0} is in C
x = x_store;
Xnew = X_fcn(x);
c= sdpvar; v = sdpvar(n,1);
exp6 = {}; L6 ={};
paras2 = [];
constrs2 = [];

for i =1:size(cx,1)
    Phi= Xnew-0.1*Xnew*cx(i,:)'*cx(i,:)*Xnew;
    [L6{i},c_L6,~] = polynomial([x;v],2);
    exp6{i} = v'*Phi*v - L6{i}*h0_x;  
    constrs2 = [constrs2 sos(exp6{i}) sos(L6{i}) c<=10];
    paras2 = [paras2; c_L6];
end
[sol2,v,Q,res] = solvesos(constrs2,-c,ops,paras2);
max_residual = max(res)
disp(sol2.info);
if sol2.problem ==0 || sol2.problem == 4 
   % ------------------------- extract X_fcn  ------------------------
    c = value(c)
    h_fcn = @(x) c- x'/X_fcn(x)*x;
end

%% plot
figure(1);clf;
step=.02;
x1= -2:step:2;
x2= -2:step/2:2;
%z=-10:step:10;
%generate a grid with all triplets (x,y,z)
[xx1,xx2] = meshgrid(x1,x2);
%intersection of inequalities in a logical matrix
% I = (xx1.*xx1 + xx2.*xx2 < 1);
I = false(size(xx1));
for i = 1:size(xx1,1)
    for j=1:size(xx1,2)
        x0 = [xx1(i,j);xx2(i,j)];
        I(i,j) = h_fcn(x0)>=0;
    end
end
ellipse(1,1,0,0,0,'r');hold on;

%plot of the points (x1,x2) that verify all inequalities
scatter(xx1(I),xx2(I),'.');
xlim([-2 2]);
axis equal

%% simulate a trajectory
x0 = [-0.6;-0.15];
dt = 0.005;
duration = 10;
% simulate
t_vec = 0:dt:duration;
T_steps = length(t_vec)-1;

% note that 
% x: theta, alpha, q
xTraj = zeros(n,T_steps+1);
uTraj = zeros(nu,T_steps);


% initial condition
x = x0;
xTraj(:,1) = x;

for i=1:T_steps
    t = t_vec(i);   

    % get the nominal state and input
    u = ustar + Y_fcn(x)/X_fcn(x)*Zx_fcn(x);

    % record
    uTraj(:,i) = u;        
    % propagate with zero-hold for control inputs
    [d_t,d_state] = ode23(@(t,x) dynamics(t,x,u),[t_vec(i) t_vec(i+1)],x); %,ode_options
    
    % update and record;
    x = d_state(end,:)';    
    xTraj(:,i+1) = x;    
end  

% plot the result
figure(1)
plot(xTraj(1,:),xTraj(2,:),'g');
xlabel('x_1');
xlabel('x_2');
figure(2);
subplot(2,1,1);
plot(t_vec,xTraj(1,:),t_vec,xTraj(2,:));
legend('x_1','x_2');
ylabel('x')
subplot(2,1,2);
plot(t_vec(1:end-1),uTraj);
ylabel('u')
%% further increase the region using the method in Li Wang's paper



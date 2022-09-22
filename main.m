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

% a standard nonlinear formulation
fx = A*x; Bx = B;   
dynamics = @(t,x,u) A*x+Bx*u;

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
deg_X = 0;
deg_Y = 4;
X_states_index = 2;

% state constraints |cx*Zx}|<1; a set C:={x||cx*Zx}|<1};
cx = x';
gx = x'*x-1; % a more general form to express safety constraints with gx<=0

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
Zx = x; Zx_store = Zx;
Zx_fcn = @(x) x;
v = sdpvar(n,1); w = sdpvar(1);
v1 = [v;w]; 
constrs = [];

% The following constraint is redundant after having constraint (5)
% (1): positive definiteness of X in D
% epsilon = 1e-5;
% [L_X,c_L_X,v_L_X] = polynomial([x;v],2);
% exp1 = v'*(X-epsilon*eye(n))*v - L_X*h0_x;
% constrs = [sos(exp1); sos(L_X)];
% paras = [paras; c_L_X];

% (1): \dot h >=0 in D
X_dot = zeros(n);
for j=1:length(X_states_index)
    X_dot = X_dot + dX_dx{j}*(A(X_states_index(j),:)*Zx);    
end
tmp = M*(A*X+B*Y); 
[Lcbf,c_Lcbf,v_Lcbf] = polynomial([x;v],2);
exp1 = v'*(X_dot-(tmp+tmp'))*v - Lcbf*h0_x;
constrs = [constrs sos(Lcbf) sos(exp1)];
paras = [paras; c_Lcbf];

% (2): |Au*u_til|<=1;
if include_input_limits 
    exp2 = {}; Lu ={};
    for j = 1:size(Au,1)    
        Psi_j = [1  Au(j,:)*Y; (Au(j,:)*Y)' X];
        [Lu{j},c_Lu,v_Lu] = polynomial([x;v1],2);
        exp2{j} = v1'*Psi_j*v1 - Lu{j}*h0_x;   
        constrs = [constrs sos(exp2{j}) sos(Lu{j})];
        paras = [paras; c_Lu];
    end
end

% (3): Ch = {x|h(x)>=0} is in C:={x||cx*Zx}|<1};
exp3 = {}; Lx ={};
for i =1:size(cx,1)
    Phi_i = [1 cx(i,:)*X; (cx(i,:)*X)' X];
    [Lx{i},c_Lx,v_Lx] = polynomial([x;v1],2);
    exp3{i} = v1'*Phi_i*v1 - Lx{i}*h0_x;  
    constrs = [constrs sos(exp3{i}) sos(Lx{i})];
    paras = [paras; c_Lx];
end

% (4): X >= X0; 
X0 = sdpvar(n);
[L_X0,c_L_X0,~] = polynomial([x;v],2);
% obj = -logdet(X0);                        % for SDPT3: no need to inlcude X0>=0 since this automatically assumed when using logdet function, https://yalmip.github.io/command/logdet/     
obj = -geomean(X0);                         % for other solvers
exp4 = v'*(X-X0)*v-L_X0*h0_x;
% constrs = [constrs  sos(L_X0) X0<=1e-3*eye(n)];             % X0>=1e-6*eye(n)
constrs = [constrs sos(exp4) sos(L_X0)];      % X0>=1e-6*eye(n)
paras = [paras; c_L_X0];

ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6); %, ,'sos.numblkdg',1e-7 sometimes this will cause infinite loop
[sol,v,Q,res] = solvesos(constrs,obj,ops,paras);
max_residual = max(res)
disp(sol.info);
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
h_fcn2 = @(x1,x2) h_fcn([x1;x2]);
%% 
figure(1);clf;
% step=.1;
% x1= -2:step:2;
% x2= -2:step/2:2;
% %generate a grid with all triplets (x,y,z)
% [xx1,xx2] = meshgrid(x1,x2);
% %intersection of inequalities in a logical matrix
% I = false(size(xx1));
% for i = 1:size(xx1,1)
%     for j=1:size(xx1,2)
%         x0 = [xx1(i,j);xx2(i,j)];
%         I(i,j) = h_fcn(x0)>=0;
%     end
% end
ellipse(1,1,0,0,0,'r');hold on;
fimplicit(h_fcn2)
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
xlim([-2 2]);
axis equal

return;

%% refine the CBF function 
% maximize c in h = c-Zx'*X^-1*Zx such that {x|h>=0} is in C
% x = x_store;
% Xnew = X_fcn(x);
% c= sdpvar; v = sdpvar(n,1);
% exp6 = {}; L6 ={};
% paras2 = [];
% constrs2 = [];
% 
% for i =1:size(cx,1)
%     Phi= Xnew-0.1*Xnew*cx(i,:)'*cx(i,:)*Xnew;
%     [L6{i},c_L6,~] = polynomial([x;v],2);
%     exp6{i} = v'*Phi*v - L6{i}*h0_x;  
%     constrs2 = [constrs2 sos(exp6{i}) sos(L6{i}) c<=10];
%     paras2 = [paras2; c_L6];
% end
% [sol2,v,Q,res] = solvesos(constrs2,-c,ops,paras2);
% max_residual = max(res)
% disp(sol2.info);
% if sol2.problem ==0 || sol2.problem == 4 
%    % ------------------------- extract X_fcn  ------------------------
%     c = value(c)
%     h_fcn = @(x) c- x'/X_fcn(x)*x;
% end

%% iteratively improve h(x)
clear v v1 L_X Lcbf Lu Lx L_X0
clear exp1 exp2 exp3 exp4 exp5
x = x_store; 
v = sdpvar(2,1);

if size(X_coef,3)> 1
    error('A constant X matrix has to be enforced so that h(x) is a polynominal function of x');    
end
% get initial hx;
P0 = eye(n)/X_coef(:,:,1);
K0 = Y_coef(:,:,1)*P0;
u_fcn0 = @(x) K0*Zx_fcn(x);
hx = 1-Zx_store'*P0*Zx_store; 

max_iter = 50;
tol = 1e-3;
min_eig_Qhs = nan*ones(1,max_iter);
epsilons = nan*ones(1,max_iter);
for k = 1:max_iter
    %% With h(x) and alpha fixed, search for u(x), Lcbf(x), and Lu_j(x)
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-7 sometimes this will cause infinite loop
    x = x_store; 
%     alpha = sdpvar; constrs =[alpha>=0 alpha<=1e2]; paras=[alpha];
    alpha = 1;  constrs =[]; paras=[];
    dh_dx = jacobian(hx,x);
    epsilon = sdpvar;
    [ux,c_u,v_u] = polynomial(x,4); % maximum: 4  best value: 4 (w/o input limits) and 

    % (1) cbf condition
    [Lcbf,c_Lcbf,v_Lcbf] = polynomial(x,2);
    exp1 = dh_dx*(fx+Bx*ux)+alpha*hx-Lcbf*hx-epsilon;
    constrs = [constrs;sos(exp1) sos(Lcbf)];
    paras = [paras;epsilon;c_u; c_Lcbf];

    % (2) control limits: |Au*u|<=1;
    if include_input_limits 
        exp1 = {}; Lu ={}; c_Lu={}; v_Lu={};
        for j = 1:size(Au,1)    
            [Lu{j},c_Lu{j},v_Lu{j}] = polynomial([x;v],2);
            exp1{j} = v'*[1  Au(j,:)*ux; (Au(j,:)*ux)' 1]*v - Lu{j}*hx;   
            constrs = [constrs sos(exp1{j}) sos(Lu{j})];
            paras = [paras; c_Lu{j}];
        end
    end
    [sol,~,Q,res] = solvesos(constrs,-epsilon,ops,paras);
    max_residual =  max(res);
    fprintf(1,'max_residual for searching ux: %.4e\n', max_residual);
    disp(sol.info)
    if max_residual > 1e-3
        disp('SOS problem not solved properly!');
%         pause;
    end
   
    if sol.problem ==0 || sol.problem == 4 
        c_u = value(c_u);
        c_Lcbf = value(c_Lcbf);
        ux = c_u'*v_u;
        Lcbf = c_Lcbf'*v_Lcbf;
        if include_input_limits 
            for j = 1:size(Au,1)  
                c_Lu{j} = value(c_Lu{j});
                Lu{j} = c_Lu{j}'*v_Lu{j};
            end
        end
    end
    epsilons(k) = value(epsilon);
    fprintf(1,'epsilon: %.3e\n',epsilons(k));

    %% With  u(x), Lcbf(x), and Lu_j(x) fixed, search for h(x)
    alpha = value(alpha);
    clear Zx
    constrs =[]; paras=[];
    Zx = monolist(x,2); % best value: 2 (w/o input limits) and 
    n_Zx = length(Zx);
    Qh = sdpvar(n_Zx);
    hx = Zx'*Qh*Zx;
    mu = sdpvar;
    dh_dx = jacobian(hx,x);
    paras = [Qh(:);mu];
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-9 sometimes this will cause infinite loop
    % (1) cbf condition
    exp1 = dh_dx*(fx+Bx*ux)+alpha*hx-Lcbf*hx;
    constrs = [Qh(1,1)>=1;Qh>=mu*eye(n_Zx);sos(exp1)];

    % (2) control limits: |Au*u|<=1;
    if include_input_limits 
        exp1 = {};
        for j = 1:size(Au,1)    
%             exp2{j} = 1-norm(Au(j,:)*ux,2)^2 - Lu{j}*hx;   
            exp1{j} = v'*[1  Au(j,:)*ux; (Au(j,:)*ux)' 1]*v - Lu{j}*hx;   
            constrs = [constrs sos(exp1{j})];
        end
    end

    % (3) state constraints: 
    for i =1:size(gx,1) 
        [Lx{i},c_Lx,v_Lx] = polynomial(x,4);
        exp2{i} = -hx-Lx{i}*gx(i);
        constrs = [constrs sos(exp2{i}) sos(Lx{i})];
        paras = [paras; c_Lx];
    end


    obj = -trace(Qh); % could be unbounded because Qh can be scaled 
    obj = -mu;
    % % obj = [];
    [sol,~,Q,res] = solvesos(constrs,obj,ops,paras);
    max_residual =  max(res);
    fprintf(1,'max_residual for searching hx: %.4e\n', max_residual);
    disp(sol.info)
    if sol.problem ==0 || sol.problem == 4 
        Qh0 = value(Qh);
%         trace_Qh = trace(Qh);
       
    end
    hx = Zx'*Qh0*Zx;  % Qh0 may be from last iteration when the iteration is failed
    min_eig_Qhs(k) = min(eig(value(Qh)));
    fprintf(1,'Minimum eigenvalue of Qh: %5.2f\n',min_eig_Qhs(k))
    if max_residual > 1e-3 ||  sol.problem ~=0 
        disp('SOS problem not solved properly. Iteration terminated!');
        break
    end
    if k>=3 && abs(min_eig_Qhs(k) - min_eig_Qhs(k-1))<tol &&  abs(min_eig_Qhs(k-1) - min_eig_Qhs(k-2))<tol
        disp('Minimum eigenvalue of Qh does not increase anymore. Iteration terminated.');     
        break;
    end
    if k== max_iter
         disp('Maximum iteration reached!');
    end
end
min_eig_Qhs = min_eig_Qhs(~isnan(min_eig_Qhs))

%% get u(x) and h(x) expressions
XX = x_store; % so that the result from "sdisplay" command uses "XX"
hx = clean(hx, 1e-10);
ux = clean(ux, 1e-10);
s_u = sdisplay(ux);
s_h = sdisplay(hx);
syms XX [n 1]
syms u_fcn [nu 1]
syms h_fcn [1 1]
for i=1:nu
    u_fcn(i,1) = eval(s_u{i});
end
h_fcn = eval(s_h{1});
% matlabFunction(u_fcn,'File','u_fcn1','Vars',{XX});
u_fcn = matlabFunction(u_fcn,'Vars',{XX});
h_fcn = matlabFunction(h_fcn,'Vars',{XX});
h_fcn2 = @(x1,x2) h_fcn([x1;x2]);

% plot
figure(1);
% step=.1;
% x1= -2:step:2;
% x2= -2:step/2:2;
% %z=-10:step:10;
% %generate a grid with all triplets (x,y,z)
% [xx1,xx2] = meshgrid(x1,x2);
% %intersection of inequalities in a logical matrix
% % I = (xx1.*xx1 + xx2.*xx2 < 1);
% I = false(size(xx1));
% for i = 1:size(xx1,1)
%     for j=1:size(xx1,2)
%         x0 = [xx1(i,j);xx2(i,j)];
%         I(i,j) = h_fcn(x0)>=0;
%     end
% end

% ellipse(1,1,0,0,0,'r');hold on;
fimplicit(h_fcn2,'c')
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
xlim([-2 2]);
axis equal
legend({'g(x)=0','h(x)=0','h_{ref}(x)=0'},'interpreter','latex');

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



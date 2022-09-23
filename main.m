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
fx = A*x; 
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

%% formulate the constraitns
Zx = x; Zx_store = Zx;
Zx_fcn = @(x) x;

plant.n = n; plant.nu = nu;
plant.A = A; plant.B = B; plant.M = M;
plant.x = x; plant.Zx = Zx;
plant.h0_x = h0_x;
plant.Au =Au; plant.cx = cx;

cbf_config.deg_X = deg_X; cbf_config.deg_Y = deg_Y;
cbf_config.X_state_index = X_states_index;
cbf_config.deg_Lcbf = 2; 
cbf_config.deg_Lu = 2;
cbf_config.deg_Lx = 2;
cbf_config.deg_L_X0 = 2;
cbf_config.include_input_limits = include_input_limits;
[h_fcn,h_fcn2,X_coef,Y_coef] = cbf_sos(plant,cbf_config)
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

%% refine the CBF function by increasing the constant 1 to c(>=1) (not needed)
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


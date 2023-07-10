%% Synthesis of CBFs using sum of squares optimization 
clear;
yalmip('clear')
%% settings
save_file = 0;
print_file = 0;
include_input_limits = 0;   % whether to include input limits
deg_X = 0;
deg_Y = 2;
X_states_index = 2;         % the index of states which X depends on
Y_states_index = [1 2];     % the index of states which Y depends on

%% dynamics and constraints
% an unstable linear system from  
% Clark, Verification and synthesis of control barrier functions, CDC, 2021
% A = [1 0; -1 4]; B=[1;0];
x = sdpvar(2,1); x_store = x; % for later use
A = [1 0; -1 0.5+x(2)]; B=[1;0];
n = size(A,1); nu = size(B,2);
zx = x; zx_store = zx; zx_fcn = @(x) x;

% written in a standard nonlinear input-affine form
fx = A*zx; 
dynamics = @(t,x,u) A*x+B*u;

% given a fixed point, rewrite the dynamics using the deviation states in linear-like form:
% fixed point
xstar = zeros(n,1);
ustar = zeros(nu,1);
% x_til = A(x_til)x_til + B(x_til)u_til, where
% x_til = x - xstar; util = u-ustar;

% x_store = x;
M = eye(n); 

% state constraints |Cx*zx}|<1; note that  = x for this example;
% a set C:={x||Cx*zx}|<1}; 
Cx = [x(1) x(2)];
cx = x'*x-1; % an alternative way to express safety constraints: qx<=0
qx_fcn = @(x,y) x^2+y^2-1;
% input constraints |Du*u|<1;
Du = 1;

% specif a set D={x|h0(x)>=0} that contains the original safety set C 
P0 = eye(n);
% In fact, D only needs to contain the resulting C_h = {x|h(x)>=0}; as a
% result, one can check the following one
% P0 = eye(n)/X_coef(:,:,1); [V,D]=eig(P0); D(1,1) = 1*D(1,1);D(n,n) = 0.1*D(n,n); P0 = V*D*V';
% plt_P0= fimplicit(@(x,y) 1-[x y]*(P0)*[x y]','b'); % check its shape;
% P0 =[1.0836   -0.6789;
%    -0.6789    6.5149];
h0_x = 1-x'*x;
h0_fcn = @(x1,x2) 1- [x1 x2]/P0*[x1 x2]';

%% formulate the constraitns


plant.n = n; plant.nu = nu; plant.N = length(zx);
plant.A = A; plant.B = B; plant.M = M;
plant.x = x; plant.zx = zx; plant.zx_fcn = zx_fcn;
plant.P0 = P0; plant.h0_x = h0_x;
plant.Du =Du; plant.Cx = Cx;

cbf_config.deg_X = deg_X; cbf_config.deg_Y = deg_Y;
cbf_config.X_state_index = X_states_index;
cbf_config.Y_state_index = Y_states_index;
cbf_config.deg_Lcbf = 2; 
cbf_config.deg_Lu = 2;
cbf_config.deg_Lx = 2;
cbf_config.deg_L_X0 = 2;
cbf_config.deg_L_P0 = 2;
cbf_config.include_input_limits = include_input_limits;
[h_init_fcn,h_init_fcn1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config);

%% 
figure(1);clf;hold on;
% old way for plotting (cumbersome)
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
%         I(i,j) = h_init_fcn(x0)>=0;
%     end
% end
% ellipse(1,1,0,0,0,'r');hold on;
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
% xlim([-2 2]);
% axis equal

% new way
fimplicit(qx_fcn,'r-','Linewidth',1); fimplicit(h0_fcn,'r--','Linewidth',1); 
plt_hinit = fimplicit(h_init_fcn1,'b--','Linewidth',1.5);
if save_file 
    file_name = ['exp1_ulim_' num2str(include_input_limits) '_degX_' num2str(deg_X) '_degY_' num2str(deg_Y) '.mat'];    
    save(file_name,'plant','cbf_config','h_init_fcn','h_init_fcn1','X_fcn','Y_fcn','X_coef','Y_coef');
end
return;

%% iteratively improve h(x)
clear v v1 L_X Lcbf Lu Lx L_X0
clear exp1 exp2 exp3 exp4 exp5
x = x_store; 
v = sdpvar(2,1);

if size(X_coef,3)> 1
    error('X has to be a constant matrix so that h(x) is a polynominal function');    
end
% get initial hx;
P0 = eye(n)/X_coef(:,:,1);
u_fcn0 = @(x) Y_fcn(x)*(P0*zx_fcn(x));
hx = 1-zx_store'*P0*zx_store; 

max_iter = 50;
tol = 1e-3;
deg_ux_refine = 4; deg_yx_4_hx_refine = 2;
min_eig_Qhs = nan*ones(1,max_iter);
epsilons = nan*ones(1,max_iter);
for k = 1:max_iter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % With h(x) and alpha() fixed, search for u(x), Lcbf(x), and Lu_j(x)
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-7 sometimes this will cause infinite loop
    x = x_store; 
%     alpha = sdpvar; constrs =[alpha>=0 alpha<=1e2]; paras=[alpha];
    alpha = 1;  constrs =[]; paras=[];
    dh_dx = jacobian(hx,x);
    epsilon = sdpvar;
    [ux,c_u,v_u] = polynomial(x,deg_ux_refine); % maximum: 4  best value: 4 (w/o input limits) and 

    % (1) cbf condition: h_dot +alpha(hx)-epsilon >= 0
    [Lcbf,c_Lcbf,v_Lcbf] = polynomial(x,2);
    exp1 = dh_dx*(fx+B*ux)+alpha*hx-Lcbf*hx-epsilon;
    constrs = [constrs;sos(exp1) sos(Lcbf)];
    paras = [paras;epsilon;c_u; c_Lcbf];

    % (2) control limits: |Du*u|<=1;
    if include_input_limits 
        exp1 = {}; Lu ={}; c_Lu={}; v_Lu={};
        for j = 1:size(Du,1)    
            [Lu{j},c_Lu{j},v_Lu{j}] = polynomial([x;v],2);
            exp1{j} = v'*[1  Du(j,:)*ux; (Du(j,:)*ux)' 1]*v - Lu{j}*hx;   
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
            for j = 1:size(Du,1)  
                c_Lu{j} = value(c_Lu{j});
                Lu{j} = c_Lu{j}'*v_Lu{j};
            end
        end
    end
    epsilons(k) = value(epsilon);
    fprintf(1,'epsilon: %.3e\n',epsilons(k));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % With  u(x), Lcbf(x), and Lu_j(x) fixed, search for h(x)
    
    alpha = value(alpha);
%     clear yx
    constrs =[]; paras=[];
    yx = monolist(x, deg_yx_4_hx_refine); 
    n_yx = length(yx);
    Qh = sdpvar(n_yx);
    hx = yx'*Qh*yx;
    mu = sdpvar;
    dh_dx = jacobian(hx,x);
    paras = [Qh(:);mu];
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-9 sometimes this will cause infinite loop
    % (1) cbf condition
    exp1 = dh_dx*(fx+B*ux)+alpha*hx-Lcbf*hx;
    constrs = [Qh(1,1)==1;Qh>=mu*eye(n_yx);sos(exp1)];

    % (2) control limits: |Du*u|<=1;
    if include_input_limits 
        exp1 = {};
        for j = 1:size(Du,1)    
%             exp2{j} = 1-norm(Du(j,:)*ux,2)^2 - Lu{j}*hx;   
            exp1{j} = v'*[1  Du(j,:)*ux; (Du(j,:)*ux)' 1]*v - Lu{j}*hx;   
            constrs = [constrs sos(exp1{j})];
        end
    end

%     (3) state constraints: ensuring that Ch={x|h(x)<=0} is a subset of
%     C{x|qx<0}, where qx = |Cx*zx|-1; in other words, the complement of C
%     is a subset of Ch. 
    for i =1:size(cx,1) 
        [Lx{i},c_Lx,v_Lx] = polynomial(x,4);
        exp2{i} = -hx-Lx{i}*cx(i);  % whenever qx(i)>0 (unsafe), we want h(x)<0 (outside CBF certified region)
        constrs = [constrs sos(exp2{i}) sos(Lx{i})];
        paras = [paras; c_Lx];
    end

    %obj = -trace(Qh); % not a good choice (could be unbounded because Qh
    %can be scaled)
    obj = -mu;
    [sol,~,Q,res] = solvesos(constrs,obj,ops,paras);
    max_residual =  max(res);
    fprintf(1,'max_residual for searching hx: %.4e\n', max_residual);
    disp(sol.info)
    if sol.problem ==0 || sol.problem == 4 
        Qh0 = value(Qh);
    end
    hx = yx'*Qh0*yx;  % Qh0 may be from last iteration when the iteration is failed
    dh_dx = jacobian(hx,x);
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
cbf_config.alpha = alpha;
cbf_config.deg_ux_refine = deg_ux_refine;
cbf_config.deg_yx_4_hx_refine = deg_yx_4_hx_refine;

%% get the expressions of u(x) and h(x) 
XX = x_store; % so that the result from "sdisplay" command uses "XX"
hx = clean(hx, 1e-10);
ux = clean(ux, 1e-10);
s_u = sdisplay(ux);
s_h = sdisplay(hx);
s_dh = sdisplay(dh_dx);
syms XX [n 1]
syms u_fcn [nu 1]
syms h_fcn [1 1]
syms dh_dx_fcn [1 n];
for i=1:nu
    u_fcn(i,1) = eval(s_u{i});
end
h_fcn = eval(s_h{1});
for i=1:n
    dh_dx_fcn(i) = eval(s_dh{i}); 
end

% matlabFunction(u_fcn,'File','u_fcn1','Vars',{XX});
u_fcn = matlabFunction(u_fcn,'Vars',{XX});
h_fcn = matlabFunction(h_fcn,'Vars',{XX});
h_fcn1 = @(x1,x2) h_fcn([x1;x2]);
dh_dx_fcn = matlabFunction(dh_dx_fcn,'Vars',{XX});
% result from Clark 2022 paper
h_fcn_existing = @(x1,x2) 1.1575-[x1-0.1378 x2]*[6.23 -26.7; -26.7 146.7]*[x1-0.1378 x2]';
if save_file 
    file_name = ['exp1_ulim_' num2str(include_input_limits) '_degY_' num2str(deg_Y) '_deg_ux_refine_' num2str(deg_ux_refine) '.mat'];    
    save(file_name,'plant','cbf_config','h_init_fcn','h_init_fcn1',...
        'h_fcn','h_fcn1','u_fcn','dh_dx_fcn','Qh0','X_fcn','Y_fcn','min_eig_Qhs');
end
	
% plot
figure(1);
% ellipse(1,1,0,0,0,'r');hold on;
plt_h = fimplicit(h_fcn1,'g','Linewidth',1.5);
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
% xlim([-2 2]);
plt_hinit = fimplicit(h_init_fcn1,'b--','Linewidth',1.5);
if include_input_limits == 0
    plt_h_existing = fimplicit(h_fcn_existing,'m','Linewidth',1.5);
    legend([plt_hinit plt_h plt_h_existing],{'$h^\textrm{init}(x)=0$','$h(x)=0$','$h(x)=0$ from [Clark 21]'},'interpreter','latex','numcolumns',2); %'Orientation','Horizontal'
else
    legend([plt_hinit plt_h],{'$h^\textrm{init}(x)=0$','$h(x)=0$'},'interpreter','latex','Orientation','Horizontal'); %'Orientation','Horizontal'
end
% axis equal
axis equal
xlim([-1.1 1.1]);
ylim([-1.3 1.3]);
xlabel('$x_1$','interpreter','latex');ylabel('$x_2$','interpreter','latex');
% legend({'$q(x)=0$','$h^\textrm{init}(x)=0$','$h(x)=0$' '$h(x)=0$ from [Clark 21]'},'interpreter','latex','Orientation','Horizontal');

goodplot([5 6]);
fig_name = ['exp1_ulim_' num2str(include_input_limits) '_cbf'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end		 

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
    u = ustar + Y_fcn(x)/X_fcn(x)*zx_fcn(x);
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



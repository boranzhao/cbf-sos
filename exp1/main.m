%% Convex Synthesis of CBFs under input constraints using SOS optimization 
clear;
yalmip('clear')
%% settings
save_file = 0;
print_file = 1;
linear_like_form = 2;       % {1,2} two options
include_input_limits = 1;   % whether to include input limits

% 
deg_X = 0; deg_Y = 2;
X_states_index = 2;         % the index of states which X depends on
Y_states_index = [1 2];     % the index of states which Y depends on

%% dynamics and constraints
% an unstable nonlinear system
n =  2; nu = 1;
x = sdpvar(n,1); x_store = x; % for later use

if linear_like_form == 1
    zx = x; zx_fcn = @(x) x;
    A = [1 0; -1 0.5+x(2)^2]; B=[1;0];
elseif linear_like_form == 2
    zx = [x;x(2)^3]; zx_fcn = @(x) [x;x(2)^3];
    A = [1 0 0; -1 0.5 1]; B=[1;0];
end
zx_store = zx; nzx = length(zx);

% written in a standard nonlinear input-affine form
fx = A*zx; 
f_fcn = @(x) [x(1); -x(1)+0.5*x(2)+x(2)^3];
dynamics = @(t,x,u) f_fcn(x)+[u;0];


% M is a matrix of N by n
M = jacobian(zx,x);  % an identity matrix when zx = x;
 % state constraint set X = {x|cx<=0}
cx = [x(1)^2-1;x(2)^2-1];
c1_fcn = @(x1,x2) x1^2-1;
c2_fcn = @(x1,x2) x2^2-1;

% (inner approximation of) state constraint set represented using zx: {x||Cx*zx|<=1};
if linear_like_form == 1 
    Cx = [x(1) 0; 0 x(2)];
elseif linear_like_form == 2
    Cx = [x(1) 0 0; 0 x(2) 0];
end

% input constraints in the form of Du0*u<=1
Du0 = 1; 

% (inner approximation of) input constraint set represented using {u||Du*u|<1};
Du = 1;

% set X_bar={x|h0(x)>=0} that contains the original safety set X
if linear_like_form == 1 
    P0 = 2*eye(nzx);
    h0_x = 1-(zx'*zx)/P0(1,1);
    h0_fcn = @(x1,x2) 1.-(x1.^2+x2.^2)/P0(1,1);
elseif linear_like_form == 2
    P0 = 3*eye(nzx);
    h0_x = 1-(zx'*zx)/P0(1,1);
    h0_fcn = @(x1,x2) 1.-(x1.^2+x2.^2+x2.^6)/P0(1,1);
    %cx = -h0_x; % uncomment this if one wants to enforce conditions in
    %X_bar instead of X
end

figure(1);clf;hold on;
plt_h0 = fimplicit(h0_fcn,'k');
fimplicit(c2_fcn,'r',[-1 1])
fimplicit(c1_fcn,'r',[-1 1])
axis equal
xlabel('$x_1$','interpreter','latex');ylabel('$x_2$','interpreter','latex')
hold on;

%% Step 1：initial CBF synthesis using the linear-like form
plant.n = n; plant.nu = nu; plant.nzx = nzx;
plant.A = A; plant.B = B; plant.M = M; plant.fx = fx; plant.dynamics=dynamics;
plant.x = x; plant.x_store = x;
plant.zx = zx; plant.zx_fcn = zx_fcn;
plant.P0 = P0; plant.h0_x = h0_x; plant.h0_fcn = h0_fcn;
plant.Du =Du; plant.Du0 = Du0;
plant.Cx = Cx; plant.cx = cx;

cbf_config.deg_X = deg_X; cbf_config.deg_Y = deg_Y;
cbf_config.X_state_index = X_states_index;
cbf_config.Y_state_index = Y_states_index;
cbf_config.deg_Lcbf = 4; 
cbf_config.deg_Lu = 2;
cbf_config.deg_Lx = 2;
cbf_config.deg_L_X0 = 2;
cbf_config.deg_L_P0 = 2; 
cbf_config.alpha = 1;
cbf_config.include_input_limits = include_input_limits;
[h_init_fcn,h_init_fcn1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config);

% plotting
figure(1);
plt_hinit = fimplicit(h_init_fcn1,'b--','Linewidth',1.5);
% saving data
if save_file 
    file_name = ['exp1_ulim_' num2str(include_input_limits) '_linearlike_' num2str(linear_like_form) '_degX_' num2str(deg_X) '_degY_' num2str(deg_Y) '.mat'];    
    save(file_name,'plant','cbf_config','h_init_fcn','h_init_fcn1','X_fcn','Y_fcn','X_coef','Y_coef');
end
return;

%% Step 2: iteratively improving the initial CBF 
clear v v1 L_X Lcbf Lu Lx L_X0
clear exp1 exp2 exp3 exp4 exp5
x = x_store; 
v = sdpvar(2,1);

if size(X_coef,3)> 1
    error('X has to be a constant matrix so that h(x) is a polynominal function');    
end
% get initial hx;
P0 = eye(nzx)/X_coef(:,:,1);
u_fcn0 = @(x) Y_fcn(x)*(P0*zx_fcn(x));
hx = 1-zx_store'*P0*zx_store; 

max_iter = 50;
tol = 1e-3;
if linear_like_form == 1
    deg_ux_refine = 3; deg_yx_4_hx_refine = 1;
elseif linear_like_form == 2
    deg_ux_refine = 5; deg_yx_4_hx_refine = 3;
end


cbf_config.deg_ux_refine = deg_ux_refine;
cbf_config.deg_yx_4_hx_refine = deg_yx_4_hx_refine;


[hx,ux,dh_dx,min_eig_Qhs] =  cbf_refine(plant,cbf_config,hx,max_iter,tol); 


%%
% get the expressions of u(x) and h(x) 
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

	
% plot
figure(1);hold on;
delete(plt_h0)
% ellipse(1,1,0,0,0,'r');hold on;
plt_h = fimplicit(h_fcn1,'g','Linewidth',1.5);
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
% xlim([-2 2]);
plt_hinit = fimplicit(h_init_fcn1,'b--','Linewidth',1.5);


% axis equal
axis equal
xlim([-1.1 1.1]);
ylim([-1.3 1.3]);
xlabel('$x_1$','interpreter','latex');ylabel('$x_2$','interpreter','latex');
% legend({'$q(x)=0$','$h^\textrm{init}(x)=0$','$h(x)=0$' '$h(x)=0$ from [Clark 21]'},'interpreter','latex','Orientation','Horizontal');


%% simulation given an initial state
min_norm_control = 0; 
x0 = [-0.5;-0.5];
if linear_like_form == 2
    x0 = [-0.8;-0.7];
end
dt = 0.005;
duration = 5;

sim_plt_save;




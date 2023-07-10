%% Convex Synthesis of CBFs under input constraints using SOS optimization 
clear;
yalmip('clear')
%% settings
save_file = 0; print_file=0;
include_input_limits = 1;

%
deg_X = 0; deg_Y = 2;
X_state_index = 1;
Y_state_index = [1 2 3];

%% dynamics and constraints 
% dynamics borroed from the following paper
% Wang, L., Han, D. and Egerstedt, M., 2018, Permissive barrier
% certificates for safe stabilization using sum-of-squares. ACC
n = 3; nu = 2;
x = sdpvar(n,1); x_store = x;            % for later use 
% linear-like form
A = [0 1 -x(3); -x(1) 0 1; -1 -2+x(2)^2 -1];
B = [0 0;1 0;0 1];
zx = x; zx_store = zx; nzx = length(zx);
zx_fcn = @(x) x;
x1=x(1); x2=x(2); x3=x(3);

% written in a standard nonlinear input-affine form
fx = [x(2)-x(3)^2;x(3)-x(1)^2;-x(1)-2*x(2)-x(3)+x(2)^3];
f_fcn = @(x) [x(2)-x(3)^2;x(3)-x(1)^2;-x(1)-2*x(2)-x(3)+x(2)^3];
dynamics = @(t,x,u) f_fcn(x)+[0;u(1);u(2)];

% % given a fixed point, rewrite the dynamics using the deviation states in linear-like form:
% % fixed point
% xstar = zeros(n,1);
% ustar = zeros(nu,1);
% % x_til = A(x_til)x_til + B(x_til)u_til, where
% % x_til = x - xstar; util = u-ustar;

% M is a matrix of nzx by n
M = jacobian(zx,x);  % an identity matrix when zx = x;

% state constraint set X = {x|cx<=0}
cx = [-(x1^2-4*x1+x2^2-2*x2+x3^2-4*x3+8);
     -(x1^2+2*x1+x2^2+4*x2+x3^2+2*x3+5);
     -(x1^2+x2^2+x3^2-12*x3+27);
     -(x1^2+x2^2+x3^2+10*x3+16)];

% (inner approximation of) state constraint set represented using zx: {x||Cx*zx|<=1};
Cx = [-x1/8+4/8 -x2/8+2/8 -x3/8+4/8;
    -x1/5-2/5 -x2/5-4/5 -x3/5-2/5;
    -x1/27 -x2/27 -x3/27+12/27;
    -x1/16 -x2/16 -x3/16-10/16];
 
syms x1 x2 x3
qx_fcn = [-(x1^2-4*x1+x2^2-2*x2+x3^2-4*x3+8);
     -(x1^2+2*x1+x2^2+4*x2+x3^2+2*x3+5);
     -(x1^2+x2^2+x3^2-12*x3+27);
     -(x1^2+x2^2+x3^2+10*x3+16)];
q1_fcn = matlabFunction(qx_fcn(1),'Vars',{x1,x2,x3});
q2_fcn = matlabFunction(qx_fcn(2),'Vars',{x1,x2,x3});
q3_fcn = matlabFunction(qx_fcn(3),'Vars',{x1,x2,x3});
q4_fcn = matlabFunction(qx_fcn(4),'Vars',{x1,x2,x3});

% input constraints in the form of Du0*u<=1
Du0 = [1 0;0 1];

% (inner approximation of) input constraint set represented using {u||Du*u|<1};
Du = [1 0;0 1];

% set X_bar={x|h0(x)>=0} that contains the original safety set X
P0 = 10*eye(3);
% h0_x = 1-0.1*(x'*x);
h0_x = 1-zx'/P0*zx;


% visualization 
figure(1);clf;
h0_fcn = @(x1,x2,x3) 1-0.01*([x1,x2,x3]*[x1,x2,x3]');
figure(1);clf;hold on;
fimplicit3(h0_fcn,'r');
fimplicit3(q1_fcn,'m','EdgeColor','k','FaceAlpha',0.5,'MeshDensity',40)
fimplicit3(q2_fcn,'FaceColor',[0.4940 0.1840 0.5560],'EdgeColor','k','FaceAlpha',0.5,'MeshDensity',40)
fimplicit3(q3_fcn,'r','EdgeColor','k','FaceAlpha',0.5)
fimplicit3(q4_fcn,'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','k','FaceAlpha',0.5)

% view(3)
view(60,30)
% xlim([-5 5]);ylim([-5 5]);zlim([-5 5]);
xlabel('$x_1$','interpreter','latex');ylabel('$x_2$','interpreter','latex');zlabel('$x_3$','interpreter','latex');

%% Step 1：initial CBF synthesis using the linear-like form
plant.n = n; plant.nu = nu; plant.nzx = nzx;
plant.A = A; plant.B = B; plant.M = M; plant.fx = fx;
plant.x = x; plant.x_store = x_store;
plant.zx = zx; plant.zx_fcn = zx_fcn;
plant.P0 = P0; plant.h0_x = h0_x; plant.cx= cx;
plant.Du =Du; plant.Du0 = Du0; plant.Cx = Cx;

cbf_config.deg_X = deg_X; cbf_config.deg_Y = deg_Y;
cbf_config.X_state_index = X_state_index;
cbf_config.Y_state_index = Y_state_index;
cbf_config.deg_Lcbf = 2; 
cbf_config.deg_Lu = 2;
cbf_config.deg_Lx = 2;
cbf_config.deg_L_X0 = 2;
cbf_config.deg_L_P0 = 2;
cbf_config.alpha = 1;
cbf_config.include_input_limits = include_input_limits;
cbf_config.q1_fcn = q1_fcn; cbf_config.q2_fcn = q2_fcn; 
cbf_config.q3_fcn = q3_fcn; cbf_config.q4_fcn = q4_fcn; 

%%%% uncomment the following line if one wants to enforce conditions in
% X_bar instead of X
plant.cx = -h0_x; 

[h_fcn_init,h_fcn_init1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config)

% plotting
figure(1);hold on;
plt_h = fimplicit3(h_fcn_init1,'facecolor',[0.3010 0.7450 0.9330],'MeshDensity',40); % close to blue [0 0.4470 0.7410] 
% saving data
if save_file 
    file_name = ['exp2_ulim_' num2str(include_input_limits) '_degX_' num2str(deg_X) '_degY_' num2str(deg_Y) '.mat'];    
    save(file_name,'plant','cbf_config','h_fcn_init','h_fcn_init1','X_fcn','Y_fcn');
end
return;

%% Step 2: iteratively improving the initial CBF 
plant.cx = cx;
clear v v1 L_X Lcbf Lu Lx L_X0
clear exp1 exp2 exp3 exp4 exp5
x = x_store; 
v = sdpvar(2,1);

if size(X_coef,3)> 1
    error('X needs to be a constant matrix so that h(x) is a polynominal function of x');    
end
% get initial hx;
P0 = eye(n)/X_coef(:,:,1);
u_fcn0 = @(x) Y_fcn(x)*(P0*zx_fcn(x));
hx = 1-zx_store'*P0*zx_store; 

max_iter = 50;
deg_ux_refine = 4; deg_yx_4_hx_refine = 2;
tol = 1e-3;
cbf_config.deg_ux_refine = deg_ux_refine;
cbf_config.deg_yx_4_hx_refine = deg_yx_4_hx_refine;

[hx,ux,dh_dx,min_eig_Qhs] =  cbf_refine(plant,cbf_config,hx,max_iter,tol); 

% get expressions of u(x) and h(x) 
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
u_fcn = matlabFunction(u_fcn,'Vars',{XX});
h_fcn = matlabFunction(h_fcn,'Vars',{XX});
h_fcn_4_plot = @(x1,x2,x3) h_fcn([x1;x2;x3]);
dh_dx_fcn = matlabFunction(dh_dx_fcn,'Vars',{XX});

if save_file 
    file_name = ['exp2_ulim_' num2str(include_input_limits) '_degY_' num2str(deg_Y) '_deg_ux_refine_' num2str(deg_ux_refine) '.mat'];    
    save(file_name,'plant','cbf_config','h_fcn_init','h_fcn_init1',...
        'h_fcn','h_fcn1','u_fcn','dh_dx_fcn','Qh0','X_fcn','Y_fcn','min_eig_Qhs');
end

% plot
figure(1);
plt_h_refine= fimplicit3(h_fcn_4_plot,'g','FaceAlpha',0.2,'LineStyle',':','EdgeColor','interp','MeshDensity',40) %'EdgeColor','interp'
% h_existing = @(x1,x2,x3) 114.3555+1.4686*x1+7.2121*x2 + 19.8479*x3-24.5412*x3^2-14.7734*x1^2-26.0129*x1*x2...
% -15.5440*x1*x3-28.3492*x2^2-27.5651*x2*x3; % from [Wang 2018, Permissive]
% h_existing = @(x1,x2,x3) 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2-13.0124; % V(x)-c from [Wang 2018, Permissive]
% plt_h_existing= fimplicit3(h_existing,'m','FaceAlpha',0.2,'LineStyle',':','EdgeColor','interp') %'EdgeColor','interp'

% axis equal
legend({'$q_1(x)=0$','$q_2(x)=0$','$q_3(x)=0$','$q_4(x)=0$','$h^\textrm{init}(x)=0$','$h(x)=0$'},'interpreter','latex');

xlim([-4 4]); ylim([-4 4]); zlim([-5 5])
view(-148.67, -13.75);

goodplot([6 8]);
set(gca, 'box', 'off');
fig_name = ['exp2_ulim_' num2str(include_input_limits) '_cbf'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end

%% simulate trajectories given initial points
close all;
min_norm_control = 0;
% options = optimoptions('quadprog','Display','off','MaxIterations',200,'LinearSolver','dense');
options = mskoptimset;

x0_range = [-3 3;-3 3; -3 3]
sim_plt_save;
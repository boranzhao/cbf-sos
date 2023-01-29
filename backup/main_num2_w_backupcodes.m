%% synthesis of CBF using SOS and linear-like form
clear;
yalmip('clear')
%%
save_file = 0; print_file=0;
% dynamics from the following paper
% Wang, L., Han, D. and Egerstedt, M., 2018, Permissive barrier
% certificates for safe stabilization using sum-of-squares. ACC
n = 3; nu = 2;
x = sdpvar(n,1); x_store = x;            % for later use 
% linear-like form
A = [0 1 -x(3); -x(1) 0 1; -1 -2+x(2)^2 -1];
B = [0 0;1 0;0 1];
% A = [0  1;-1 0]; B = [0 1];
% x = sdpvar(2,1);
zx = x; zx_store = zx;
zx_fcn = @(x) x;
x1=x(1); x2=x(2); x3=x(3);


% a standard nonlinear form
fx = [x(2)-x(3)^2;x(3)-x(1)^2;-x(1)-2*x(2)-x(3)+x(2)^3];
f_fcn = @(x) [x(2)-x(3)^2;x(3)-x(1)^2;-x(1)-2*x(2)-x(3)+x(2)^3];
dynamics = @(t,x,u) f_fcn(x)+[0;u(1);u(2)];

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
deg_Y = 2;
X_states_index = 1;

% state constraints |cx*zx}|<1; a set C:={x||cx*zx}|<1};
cx = [-x1/8+4/8 -x2/8+2/8 -x3/8+4/8;
    -x1/5-2/5 -x2/5-4/5 -x3/5-2/5;
    -x1/27 -x2/27 -x3/27+12/27;
    -x1/16 -x2/16 -x3/16-10/16];
 % a more general form to express safety constraints with qx<=0
qx = [-(x1^2-4*x1+x2^2-2*x2+x3^2-4*x3+8);
     -(x1^2+2*x1+x2^2+4*x2+x3^2+2*x3+5);
     -(x1^2+x2^2+x3^2-12*x3+27);
     -(x1^2+x2^2+x3^2+10*x3+16)];
syms x1 x2 x3
qx_fcn = [-(x1^2-4*x1+x2^2-2*x2+x3^2-4*x3+8);
     -(x1^2+2*x1+x2^2+4*x2+x3^2+2*x3+5);
     -(x1^2+x2^2+x3^2-12*x3+27);
     -(x1^2+x2^2+x3^2+10*x3+16)];
q1_fcn = matlabFunction(qx_fcn(1),'Vars',{x1,x2,x3});
q2_fcn = matlabFunction(qx_fcn(2),'Vars',{x1,x2,x3});
q3_fcn = matlabFunction(qx_fcn(3),'Vars',{x1,x2,x3});
q4_fcn = matlabFunction(qx_fcn(4),'Vars',{x1,x2,x3});
% input constraints |Au*u|<1;
Au = [1 0;0 1];

% specif a set D={x|h0(x)>=0} that contains the original safety set
h0_x = 1-0.1*(x'*x);

%% visualization 
figure(1);clf;
h0_fcn = @(x1,x2,x3) 1-0.01*([x1,x2,x3]*[x1,x2,x3]');
figure(1);clf;hold on;
fimplicit3(q1_fcn,'m','EdgeColor','k','FaceAlpha',0.5,'MeshDensity',40)
fimplicit3(q2_fcn,'FaceColor',[0.4940 0.1840 0.5560],'EdgeColor','k','FaceAlpha',0.5,'MeshDensity',40)
fimplicit3(q3_fcn,'r','EdgeColor','k','FaceAlpha',0.5)
fimplicit3(q4_fcn,'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor','k','FaceAlpha',0.5)
% fimplicit3(h0_fcn,'r');
% axis equal
% view(3)
view(60,30)
% xlim([-5 5]);ylim([-5 5]);zlim([-5 5]);
xlabel('x_1','interpreter','latex');ylabel('x_2','interpreter','latex');zlabel('x_3','interpreter','latex');

%% formulate the constraitns
plant.n = n; plant.nu = nu;
plant.A = A; plant.B = B; plant.M = M;
plant.x = x; plant.zx = zx;
plant.h0_x = h0_x;
plant.Au =Au; plant.cx = cx;

cbf_config.deg_X = deg_X; cbf_config.deg_Y = deg_Y;
cbf_config.X_state_index = X_states_index;
cbf_config.deg_Lcbf = 2; 
cbf_config.deg_Lu = 2;
cbf_config.deg_Lx = 2;
cbf_config.deg_L_X0 = 2;
cbf_config.include_input_limits = include_input_limits;
cbf_config.q1_fcn = q1_fcn; cbf_config.q2_fcn = q2_fcn; 
cbf_config.q3_fcn = q3_fcn; cbf_config.q4_fcn = q4_fcn; 
[h_fcn_init,h_fcn_init1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config)
%% 
figure(1);
% old way
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
%         I(i,j) = h_fcn_init(x0)>=0;
%     end
% end
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
% axis equal

% new way
plt_h = fimplicit3(h_fcn_init1,'facecolor',[0.3010 0.7450 0.9330],'MeshDensity',40); % close to blue [0 0.4470 0.7410]

    
if save_file 
    file_name = ['exp2_ulim_' num2str(include_input_limits) '_degX_' num2str(deg_X) '_degY_' num2str(deg_Y) '.mat'];    
    save(file_name,'plant','cbf_config','h_fcn_init','h_fcn_init1','X_fcn','Y_fcn');
end
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
u_fcn0 = @(x) K0*zx_fcn(x);
hx = 1-zx_store'*P0*zx_store; 

max_iter = 50;
deg_ux_refine = 4; deg_yx_4_hx_refine = 2;
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
   
    [ux1,c_u1,v_u] = polynomial(x,deg_ux_refine); %
    c_u2 = sdpvar(length(c_u1),1); ux2 = c_u2'*v_u;
    ux = [ux1;ux2];c_u=[c_u1;c_u2];
    
    % (1) cbf condition: h_dot +alpha(hx)-epsilon >= 0
    [Lcbf,c_Lcbf,v_Lcbf] = polynomial(x,2);
    exp1 = dh_dx*(fx+B*ux)+alpha*hx-Lcbf*hx-epsilon;
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
        c_u1 = value(c_u1); c_u2 = value(c_u2);
        c_Lcbf = value(c_Lcbf);
        ux = [c_u1 c_u2]'*v_u;
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
    constrs =[]; paras=[];
    yx = monolist(x,deg_yx_4_hx_refine); % best value: 2 (w/o input limits) and 
    n_yx = length(yx);
    Qh = sdpvar(n_yx);
    hx = yx'*Qh*yx;
    mu = sdpvar;
    dh_dx = jacobian(hx,x);
    paras = [Qh(:);mu];
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-9 sometimes this will cause infinite loop
    % (1) cbf condition
    exp1 = dh_dx*(fx+B*ux)+alpha*hx-Lcbf*hx;
    constrs = [Qh(1,1)>=1;Qh>=mu*eye(n_yxh);sos(exp1)];

    % (2) control limits: |Au*u|<=1;
    if include_input_limits 
        exp1 = {};
        for j = 1:size(Au,1)    
%             exp2{j} = 1-norm(Au(j,:)*ux,2)^2 - Lu{j}*hx;   
            exp1{j} = v'*[1  Au(j,:)*ux; (Au(j,:)*ux)' 1]*v - Lu{j}*hx;   
            constrs = [constrs sos(exp1{j})];
        end
    end

    % (3) state constraints: ensuring that Ch={x|h(x)<=0} is a subset of
%     C{x|qx<0}, where qx = |cx*zx|-1; in other words, the complement of C
%     is a subset of Ch. 
    for i =1:size(qx,1) 
        [Lx{i},c_Lx,v_Lx] = polynomial(x,4);
        exp2{i} = -hx-Lx{i}*qx(i);  %whenever qx(i)>0 (unsafe), we want h(x)<0 (outside CBF certified region)
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

%% get expressions of u(x) and h(x) 
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
h_fcn_refine1 = @(x1,x2,x3) h_fcn([x1;x2;x3]);
dh_dx_fcn = matlabFunction(dh_dx_fcn,'Vars',{XX});

if save_file 
    file_name = ['exp2_ulim_' num2str(include_input_limits) '_degY_' num2str(deg_Y) '_deg_ux_refine_' num2str(deg_ux_refine) '.mat'];    
    save(file_name,'plant','cbf_config','h_init_fcn','h_init_fcn1',...
        'h_fcn','h_fcn1','u_fcn','dh_dx_fcn','Qh0','X_fcn','Y_fcn','min_eig_Qhs');
end

%% plot
figure(1);hold on;
plt_h_refine= fimplicit3(h_fcn1,'g','FaceAlpha',0.2,'LineStyle',':','EdgeColor','interp') %'EdgeColor','interp'
% h_existing = @(x1,x2,x3) 114.3555+1.4686*x1+7.2121*x2 + 19.8479*x3-24.5412*x3^2-14.7734*x1^2-26.0129*x1*x2...
% -15.5440*x1*x3-28.3492*x2^2-27.5651*x2*x3; % from [Wang 2018, Permissive]
h_existing = @(x1,x2,x3) 5*x1^2+10*x1*x2+2*x1*x3+10*x2^2+6*x2*x3+4*x3^2-13.0124; % V(x)-c from [Wang 2018, Permissive]
plt_h_existing= fimplicit3(h_existing,'m','FaceAlpha',0.2,'LineStyle',':','EdgeColor','interp') %'EdgeColor','interp'

% axis equal
legend({'$q_1(x)=0$','$q_2(x)=0$','$q_3(x)=0$','$q_4(x)=0$','$h(x)=0$','$h^\textrm{ref}(x)=0$'},'interpreter','latex');

view(-148.67, -13.75);

goodplot([6 8]);
set(gca, 'box', 'off');
fig_name = ['exp2_ulim_' num2str(include_input_limits) '_cbf'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end

%% simulate trajectories given initial points
min_norm_control = 1;
% options = optimoptions('quadprog','Display','off','MaxIterations',200,'LinearSolver','dense');
options = mskoptimset;
if include_input_limits
    step = 0.5;
else
    step=1;
end
x1= -3:step:3;
x2= -3:step:3;
x3= -3:step:3;
[xx1,xx2,xx3] = meshgrid(x1,x2,x3);
%intersection of inequalities in a logical matrix
% I = (xx1.*xx1 + xx2.*xx2 < 1);
I = false(size(xx1));
x0s = [];
for i = 1:size(xx1,1)
    for j=1:size(xx1,2)
        for k= 1:size(xx1,3)
            x0 = [xx1(i,j,k);xx2(i,j,k);xx3(i,j,k)];
            if h_fcn(x0)>=0 
                if isempty(x0s) || (~isempty(x0s) && all(sum((x0-x0s).^2)>1))
                   I(i,j,k) = 1;
                   x0s = [x0s x0];
                end
            end
        end
    end
end
% x0 = [-1;3;-1]; 
% x0 = [2;-1;0];
figure(2);clf;hold on;
plt_h_refine= fimplicit3(h_fcn1,'g','FaceAlpha',0.15,'LineStyle',':','EdgeColor','interp'); %'EdgeColor','interp'
figure(3);clf;hold on;
for i = 1:size(x0s,2)
    x0 = x0s(:,i);
    if h_fcn(x0)<0
        error('Initial state is outside the safe region');
    end
    dt = 0.005;
    duration = 20;
    % simulate
    t_vec = 0:dt:duration;
    T_steps = length(t_vec)-1;

    xTraj = zeros(n,T_steps+1);
    uTraj = zeros(nu,T_steps);

    % initial condition
    x = x0;
    xTraj(:,1) = x;

    for i=1:T_steps
        t = t_vec(i);   
        % get the nominal state and input
%         u = ustar + Y_fcn(x)/X_fcn(x)*zx_fcn(x); % original control law
        if min_norm_control == 0
            u  = u_fcn(x);                             % refined control law
        else
            % min-norm control law with the searched cbf
            dh_dx0 = dh_dx_fcn(x);
            % phi0 + phi1*u >=0 
            phi0 = dh_dx0*f_fcn(x) + alpha*h_fcn(x);
            phi1 = dh_dx0*B;  
            % without constraint
%             if phi0>=0
%                 u = zeros(nu,1);
%             else
%                 u = -phi1'*phi0/(phi1*phi1');
%             end 
            
            % additionally add constraint            
            [u,feval,exitflag,~] = quadprog(eye(nu),[0 0]',-phi1,phi0,[],[],-1./diag(Au),1./diag(Au),[],options); % options
        end
       
        % record
        uTraj(:,i) = u;        
        % propagate with zero-hold for control inputs
        [d_t,d_state] = ode23(@(t,x) dynamics(t,x,u),[t_vec(i) t_vec(i+1)],x); %,ode_options

        % update and record;
        x = d_state(end,:)';    
        xTraj(:,i+1) = x;    
    end  

    % plot the result
    figure(2)
    plot3(xTraj(1,:),xTraj(2,:),xTraj(3,:),'k','Linewidth',1.5);
    scatter3(xTraj(1,1),xTraj(2,1),xTraj(3,1),'k','Linewidth',1.5);
    scatter3(xTraj(1,end),xTraj(2,end),xTraj(3,end),'b*','Linewidth',1.5);
    figure(3);
%     subplot(2,1,1);hold on;
%     h_plt_x1 = plot(t_vec,xTraj(1,:),'k','Linewidth',1);
%     h_plt_x2 = plot(t_vec,xTraj(2,:),'b','Linewidth',1);
%     h_plt_x3 = plot(t_vec,xTraj(3,:),'r','Linewidth',1);
%     
%     subplot(2,1,2); hold on;
    h_plt_u1 = plot(t_vec(1:end-1),uTraj(1,:),'k','Linewidth',1);
    h_plt_u2 = plot(t_vec(1:end-1),uTraj(2,:),'r','Linewidth',1);
end
figure(2);
scatter3(0,0,0,200,'r','filled','p','MarkerFaceColor','r');
fig_name = ['exp2_ulim_' num2str(include_input_limits) '_min_norm_', num2str(min_norm_control) '_3d_trajs']; 
view(130.2,37.8);
xlabel('x_1','interpreter','latex');ylabel('x_2','interpreter','latex');zlabel('x_3','interpreter','latex');
goodplot([6 6]);
legend off;
set(gca, 'box', 'off');
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end
figure(3)
% subplot(2,1,1)
% legend([h_plt_x1 h_plt_x2 h_plt_x3],{'x_1','x_2','x_3'});
% ylabel('x')
% goodplot([6 6]);
% subplot(2,1,2)
ylabel('u','interpreter','latex')
xlabel('Time (s)');
legend([h_plt_u1 h_plt_u2],{'u_1','u_2'},'Orientation','Horizontal','interpreter','latex');
goodplot([6 6]);
fig_name = ['exp2_ulim_' num2str(include_input_limits) '_min_norm_', num2str(min_norm_control) '_trajs']; 
goodplot([6 3]);
% set(gca, 'box', 'off');
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end

%% refine the CBF function by increasing the constant 1 to c(>=1) (not needed)
% maximize c in h = c-zx'*X^-1*zx such that {x|h>=0} is in C
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
%     h_init_fcn = @(x) c- x'/X_fcn(x)*x;
% end


%% for plot (backup)
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
%         I(i,j) = h_init_fcn(x0)>=0;
%     end
% end
%plot of the points (x1,x2) that verify all inequalities
% scatter(xx1(I),xx2(I),'.');
% xlim([-2 2]);
function [h_fcn,h_fcn1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config)

n = plant.n; nu = plant.nu; nzx = plant.nzx;
A = plant.A; B = plant.B; M = plant.M;
x = plant.x; zx = plant.zx; zx_fcn = plant.zx_fcn;
P0 = plant.P0; %h0_x = plant.h0_x;
Du = plant.Du;  Cx = plant.Cx; cx = plant.cx; %poly_Cx = plant.poly_Cx; 

deg_X = cbf_config.deg_X; deg_Y = cbf_config.deg_Y;
X_state_index = cbf_config.X_state_index;
Y_state_index = cbf_config.Y_state_index;
deg_Lcbf = cbf_config.deg_Lcbf; 
deg_Lu = cbf_config.deg_Lu;
deg_Lx = cbf_config.deg_Lx;
deg_L_X0 = cbf_config.deg_L_X0;
deg_L_P0 = cbf_config.deg_L_P0;

% polynomials and matrix coefficients
polyX = monolist(x(X_state_index),deg_X); % note that X can only depend on the states whose derivatives are not influenced by control inputs
polyY = monolist(x(Y_state_index),deg_Y);    
d_polyX = jacobian(polyX,x(X_state_index));
X_coef = sdpvar(nzx,nzx,length(polyX),'symmetric');
Y_coef = sdpvar(nu,nzx,length(polyY));
X = zeros(nzx); Y = zeros(nu,nzx);
for j=1:length(X_state_index)
    dX_dx{j} = sdpvar(nzx,nzx);
end
for i = 1:length(polyX)
    X = X + X_coef(:,:,i)*polyX(i);
    for j = 1:length(X_state_index)
        dX_dx{j} =  dX_dx{j} - dX_dx{j}*(i==1) + X_coef(:,:,i)*d_polyX(i,j);
    end    
end

for i = 1:length(polyY)
    Y = Y + Y_coef(:,:,i)*polyY(i);
end
paras = [X_coef(:);Y_coef(:)];
v = sdpvar(nzx,1); w = sdpvar(1);
v1 = [v;w]; 
constrs = [];

% (1): \dot h >=0 for any x in X
% Q: why not use \dot h + alpha*h >= 0?  A: when we use the linear-like
% form, we can only impose \dot h>=0 to get tractable solutions. 
X_dot = zeros(nzx);
for j=1:length(X_state_index)
    X_dot = X_dot + dX_dx{j}*(A(X_state_index(j),:)*zx);    
end
tmp = M*(A*X+B*Y); 

exp1 = v'*(X_dot-(tmp+tmp'))*v;
for i = 1:length(cx)
    [Lcbf,c_Lcbf,~] = polynomial([x;v],deg_Lcbf); 
    exp1 = exp1- Lcbf*(-cx(i));
    constrs = [constrs sos(Lcbf)];
    paras = [paras; c_Lcbf];
end
constrs = [constrs sos(exp1)];

% (2): |Du*u_til|<=1 for any x in X
if cbf_config.include_input_limits 
    exp2 = {};
    for j = 1:size(Du,1)    
        Psi_j = [1  Du(j,:)*Y; (Du(j,:)*Y)' X];
        exp2{j} = v1'*Psi_j*v1;
        for i = 1:length(cx)
            [Lu,c_Lu,~] = polynomial([x;v1],deg_Lu);
            exp2{j} = exp2{j}- Lu*(-cx(i));
            constrs = [constrs sos(Lu)];
            paras = [paras; c_Lu];
        end
        constrs = [constrs sos(exp2{j})];
    end
end

% (3): Xh = {x|h(x)>=0} is a subset of X_ubar:={x||cx*zx}|<1} for any x in X
exp3 = {}; Lx ={};
for i =1:size(Cx,1)
    Phi_i = [1 Cx(i,:)*X; (Cx(i,:)*X)' X];
    exp3{i} = v1'*Phi_i*v1;  
    for k = 1:length(cx)
        [Lx,c_Lx,~] = polynomial([x;v1],deg_Lx);
        exp3{i} = exp3{i}- Lx*(-cx(k));
        constrs = [constrs  sos(Lx)];
        paras = [paras; c_Lx];
    end
    constrs = [constrs sos(exp3{i})];
end

% (4): X>=X0 for any x in X
X0 = sdpvar(nzx);
exp4 = v'*(X-X0)*v;
for k = 1:length(cx)
    [L_X0,c_L_X0,~] = polynomial([x;v],deg_L_X0); 
    exp4 = exp4 - L_X0*(-cx(k));
    constrs = [constrs  sos(L_X0)];
    paras = [paras; c_L_X0];
end
constrs = [constrs sos(exp4)];


% exp4 = v'*(X-X0)*v-L_X0*h0_x;
% constrs = [constrs sos(exp4) sos(L_X0)];   %X0>=1e-3*eye(n)     
% paras = [paras; c_L_X0]; 

% (5): P0>=X (so that Xh is a subset of X_bar) for any x in X              
exp5 = v'*(P0-X)*v;
for k = 1:length(cx)
    [L_P0,c_L_P0,~] = polynomial([x;v],deg_L_P0);
    exp5 = exp5 - L_P0*(-cx(k));
    constrs = [constrs  sos(L_P0)];
    paras = [paras; c_L_P0];
end
constrs = [constrs sos(exp5)];        

% cost function
% no need to inlcude X0>=0 when using logdet or geomean function, https://yalmip.github.io/command/logdet/ 
% using logdet for SDPT3 and geomean for other solvers such as sedumi, mosek, etc. 
% obj = -logdet(X0);            
obj = -geomean(X0);                           

ops = sdpsettings('solver','mosek','sos.numblkdg',1e-9); % sometimes need to change 'sos.numblkdg' to avoid an infinite loop
[sol,v,Q,res] = solvesos(constrs,obj,ops,paras);
max_residual = max(res)
disp(sol.info);
if sol.problem ==0 || sol.problem == 4 
   % ------------------------- extract X_fcn  ------------------------
    X_coef = value(X_coef) 
    X_coef(abs(X_coef)<=1e-10) = 0;
    
    Y_coef = value(Y_coef);
    X_fcn = zeros(nzx); Y_fcn = zeros(nu,nzx);
    for i=1:length(polyX)
        X_fcn = X_fcn+ X_coef(:,:,i)*polyX(i);
    end
    for i=1:length(polyY)
        Y_fcn = Y_fcn+ Y_coef(:,:,i)*polyY(i);
    end
    X_fcn = clean(X_fcn, 1e-10);       
    Y_fcn = clean(Y_fcn, 1e-10);
    XX = plant.x; % so that the result from "sdisplay" command uses "XX"
    s = sdisplay(X_fcn);
    s2 = sdisplay(Y_fcn);
    clear X_fcn Y_fcn
    syms XX [n 1]
    syms X_fcn [nzx nzx]
    Y_fcn = sym('Y_fcn',[nu nzx]);
    for i=1:nzx
        for j=1:nzx
            X_fcn(i,j) = eval(s{i,j});
        end
    end
    X_fcn = matlabFunction(X_fcn,'Vars',{XX});
    
    for i=1:nu
        for j=1:nzx
            Y_fcn(i,j) = eval(s2{i,j});
        end
    end
    Y_fcn = matlabFunction(Y_fcn,'Vars',{XX});
end
h_fcn = @(x) 1- zx_fcn(x)'/X_fcn(x)*zx_fcn(x);
switch length(x)
    case 2
        h_fcn1 = @(x1,x2) h_fcn([x1;x2]); 
    case 3
        h_fcn1 = @(x1,x2,x3) h_fcn([x1;x2;x3]);
    case 4
        h_fcn1 = @(x1,x2,x3,x4) h_fcn([x1;x2;x3;x4]); 
    case 5
        h_fcn1 = @(x1,x2,x3,x4,x5) h_fcn([x1;x2;x3;x4;x5]);
    case 6
        h_fcn1 = @(x1,x2,x3,x4,x5,x6) h_fcn([x1;x2;x3;x4;x5;x6]); 
end
    
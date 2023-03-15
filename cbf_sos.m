function [h_fcn,h_fcn1,X_coef,Y_coef,X_fcn,Y_fcn] = cbf_sos(plant,cbf_config)
% extracting parameters
n= plant.n; nu = plant.nu;
A = plant.A; B = plant.B; M = plant.M;
x = plant.x; zx = plant.zx;
P0 = plant.P0; h0_x = plant.h0_x;
Du = plant.Du; cx = plant.cx;

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
X_coef = sdpvar(n,n,length(polyX),'symmetric');
Y_coef = sdpvar(nu,n,length(polyY));
X = zeros(n); Y = zeros(nu,n);
% dX_dx = zeros(n,n,length(X_states_index));
for j=1:length(X_state_index)
    dX_dx{j} = sdpvar(n,n);
end
% dX_dx = X_coef(:,:,1)*d_polyX(1,:);
% dX_dx(dX_dx-dX_dx==0)=0;
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
v = sdpvar(n,1); w = sdpvar(1);
v1 = [v;w]; 
constrs = [];

% (1): \dot h >=0 in D:
% Q: why not use \dot h + alpha*h >= 0?  A: when we use the linear-like
% form, we can only impose \dot h>=0 to get tractable solutions. 
X_dot = zeros(n);
for j=1:length(X_state_index)
    X_dot = X_dot + dX_dx{j}*(A(X_state_index(j),:)*zx);    
end
tmp = M*(A*X+B*Y); 
[Lcbf,c_Lcbf,v_Lcbf] = polynomial([x;v],deg_Lcbf); %%%%%%%%%%%%%%%%%%%%%% this may need to be changed
exp1 = v'*(X_dot-(tmp+tmp'))*v - Lcbf*h0_x;
constrs = [constrs sos(Lcbf) sos(exp1)];
paras = [paras; c_Lcbf];

% (2): |Du*u_til|<=1 (input constraint) in D
if cbf_config.include_input_limits 
    exp2 = {}; Lu ={};
    for j = 1:size(Du,1)    
        Psi_j = [1  Du(j,:)*Y; (Du(j,:)*Y)' X];
        [Lu{j},c_Lu,v_Lu] = polynomial([x;v1],deg_Lu);
        exp2{j} = v1'*Psi_j*v1 - Lu{j}*h0_x;   
        constrs = [constrs sos(exp2{j}) sos(Lu{j})];
        paras = [paras; c_Lu];
    end
end

% (3): Ch = {x|h(x)>=0} is a subset of C:={x||cx*zx}|<1};
exp3 = {}; Lx ={};
for i =1:size(cx,1)
    Phi_i = [1 cx(i,:)*X; (cx(i,:)*X)' X];
    [Lx{i},c_Lx,v_Lx] = polynomial([x;v1],deg_Lx);
    exp3{i} = v1'*Phi_i*v1 - Lx{i}*h0_x;  
    constrs = [constrs sos(exp3{i}) sos(Lx{i})];
    paras = [paras; c_Lx];
end

% (4): X >= X0; 
X0 = sdpvar(n);
[L_X0,c_L_X0,~] = polynomial([x;v],deg_L_X0);                        
exp4 = v'*(X-X0)*v-L_X0*h0_x;
constrs = [constrs sos(exp4) sos(L_X0)];       
paras = [paras; c_L_X0];

% (5): P0>=X in D so that Ch is a subset of D
[L_P0,c_L_P0,~] = polynomial([x;v],deg_L_P0);                        
exp5 = v'*(P0-X)*v-L_P0*h0_x;
constrs = [constrs sos(exp5) sos(L_P0)];        
paras = [paras; c_L_P0];


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
    X_fcn = zeros(n); Y_fcn = zeros(nu,n);
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
    syms XX [n 1]
    syms X_fcn [n n]
    syms Y_fcn [nu n]
    for i=1:n
        for j=1:n
            X_fcn(i,j) = eval(s{i,j});
        end
    end
%     matlabFunction(X_fcn,'File','X_fcn1','Vars',{XX});
    X_fcn = matlabFunction(X_fcn,'Vars',{XX});
    
    for i=1:nu
        for j=1:n
            Y_fcn(i,j) = eval(s2{i,j});
        end
    end
%     matlabFunction(Y_fcn,'File','Y_fcn1','Vars',{XX});
    Y_fcn = matlabFunction(Y_fcn,'Vars',{XX});
end
h_fcn = @(x) 1- x'/X_fcn(x)*x;
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
    
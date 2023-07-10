function  [hx,ux,dh_dx,min_eig_Qhs] =  cbf_refine(plant,cbf_config,hx,max_iter,tol)
x = plant.x_store; 
%     alpha = sdpvar; constrs =[alpha>=0 alpha<=1e2]; paras=[alpha];
alpha = cbf_config.alpha;
deg_ux_refine = cbf_config.deg_ux_refine;
deg_yx_4_hx_refine = cbf_config.deg_yx_4_hx_refine;

min_eig_Qhs = nan*ones(1,max_iter);
epsilons = nan*ones(1,max_iter);
for k = 1:max_iter
    %% With h(x) and alpha fixed, search for u(x), Lcbf(x), and Lu_j(x)
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-7 sometimes this will cause infinite loop
    constrs =[]; paras=[];
    dh_dx = jacobian(hx,x);
    epsilon = sdpvar;    
   
    [ux,c_u,v_u] = polynomial(x,deg_ux_refine); %
    if plant.nu>1
        for j = 2:plant.nu
            c_uj = sdpvar(length(c_u),1); uxj = c_uj'*v_u;
            ux = [ux;uxj];c_u=[c_u c_uj];
        end
    end
    
    % (1) cbf condition: h_dot +alpha(hx)-epsilon >= 0
    [Lcbf,c_Lcbf,v_Lcbf] = polynomial(x,cbf_config.deg_Lcbf);
    exp1 = dh_dx*(plant.fx+plant.B*ux)+alpha*hx-Lcbf*hx-epsilon;
    constrs = [constrs;sos(exp1) sos(Lcbf)];
    paras = [paras;epsilon;c_u(:); c_Lcbf];

    % (2) control limits: Du0*u <=1;
    if cbf_config.include_input_limits 
        exp1 = {}; Lu ={}; c_Lu={}; v_Lu={};
        for j = 1:size(plant.Du0,1)
            [Lu{j},c_Lu{j},v_Lu{j}] = polynomial(x,cbf_config.deg_Lu); 
            exp1{j} = 1-plant.Du0(j,:)*ux- Lu{j}*hx;   % Du0(j,:)*ux<=1 for any x in Ch
            
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
    end
   
    if sol.problem ==0 || sol.problem == 4 
        c_u = value(c_u); %c_u2 = value(c_u2);
        c_Lcbf = value(c_Lcbf);
        ux = c_u'*v_u;
        Lcbf = c_Lcbf'*v_Lcbf;
        if cbf_config.include_input_limits 
            for j = 1:size(plant.Du0,1)  
                c_Lu{j} = value(c_Lu{j});
                Lu{j} = c_Lu{j}'*v_Lu{j};
            end
        end
    end
    epsilons(k) = value(epsilon);
    fprintf(1,'epsilon: %.3e\n',epsilons(k));

    %% With  u(x), Lcbf(x), and Lu_j(x) fixed, search for h(x)
    constrs =[]; paras=[];
    yx = monolist(x,deg_yx_4_hx_refine); 
    n_yx = length(yx);
    Qh = sdpvar(n_yx);
    hx = yx'*Qh*yx;
    mu = sdpvar;
    dh_dx = jacobian(hx,x);
    paras = [Qh(:);mu];
    ops = sdpsettings('solver','mosek','sos.numblkdg',1e-6,'verbose',0); %, ,'sos.numblkdg',1e-9 sometimes this will cause infinite loop
    % (1) cbf condition
    exp1 = dh_dx*(plant.fx+plant.B*ux)+alpha*hx-Lcbf*hx;
    constrs = [Qh(1,1)==1;Qh>=mu*eye(n_yx);sos(exp1)];

    % (2) control limits: Du0*u<=1;
    if cbf_config.include_input_limits 
        exp1 = {};
        for j = 1:size(plant.Du0,1)    
            exp1{j} = 1-plant.Du0(j,:)*ux- Lu{j}*hx;
            constrs = [constrs sos(exp1{j})];
        end
    end

    % (3) state constraints: ensuring that Xh={x|h(x)<=0} is a subset of
%     X={x|cx<=0}, in other words, the complement of X
%     is a subset of the complement of Xh. 
    for i =1:size(plant.cx,1) 
        [Lx{i},c_Lx,v_Lx] = polynomial(x,cbf_config.deg_Lx); %may need to manually check wheter the degree is appropriate
        exp2{i} = -hx-Lx{i}*plant.cx(i);  %whenever cx(i)>0 (unsafe), we want h(x)<0 (outside CBF certified region)
        constrs = [constrs sos(exp2{i}) sos(Lx{i})];
        paras = [paras; c_Lx];
    end

    obj = -mu;
    [sol,~,Q,res] = solvesos(constrs,obj,ops,paras);
    max_residual =  max(res);
    fprintf(1,'max_residual for searching hx: %.4e\n', max_residual);
    disp(sol.info)
    if sol.problem ==0 || sol.problem == 4 
        Qh0 = value(Qh);   
    end
    hx = yx'*Qh0*yx;  % Qh0 may be from last iteration when the iteration failed
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
min_eig_Qhs = min_eig_Qhs(~isnan(min_eig_Qhs));
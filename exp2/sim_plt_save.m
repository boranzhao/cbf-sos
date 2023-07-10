color = {'k','b','r','m','c','g',[0 0.4470 0.7410],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]}
n_color = length(color);
step = 1.2;
x1= x0_range(1,1):step:x0_range(1,2);
x2= x0_range(2,1):step:x0_range(2,2);
x3= x0_range(3,1):step:x0_range(3,2);
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
% x0 = [2;-1;0];
figure(2);clf;hold on;
plt_h_refine= fimplicit3(h_fcn_4_plot,'g','FaceAlpha',0.15,'LineStyle',':','EdgeColor','interp','MeshDensity',20) %'EdgeColor','interp'
figure(3);clf;
%hold on;
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

    for k=1:T_steps
        t = t_vec(k);   
        % get the nominal state and input
%         u = ustar + Y_fcn(x)/X_fcn(x)*zx_fcn(x); % original control law
        if min_norm_control == 0
            u  = u_fcn(x);                             % refined control law
        else
            % min-norm control law with the searched cbfh_fcn1
            dh_dx0 = dh_dx_fcn(x);
            % phi0 + phi1*u >=0 
            phi0 = dh_dx0*f_fcn(x) + cbf_config.alpha*h_fcn(x);
            phi1 = dh_dx0*B;  
            % without constraint
%             if phi0>=0
%                 u = zeros(nu,1);
%             else
%                 u = -phi1'*phi0/(phi1*phi1');
%             end 
            
            % additionally add constraint            
            [u,feval,exitflag,~] = quadprog(eye(nu),[0 0]',-phi1,phi0,[],[],-1./diag(Du),1./diag(Du),[],options); % options
        end
       
        % record
        uTraj(:,k) = u;        
        % propagate with zero-hold for control inputs
        [d_t,d_state] = ode23(@(t,x) dynamics(t,x,u),[t_vec(k) t_vec(k+1)],x); %,ode_options

        % update and record;
        x = d_state(end,:)';    
        xTraj(:,k+1) = x;    
    end  

    % plot the result
    figure(2)
    plot3(xTraj(1,:),xTraj(2,:),xTraj(3,:),'color',color{mod(i-1,n_color)+1},'Linewidth',1.5);
    scatter3(xTraj(1,1),xTraj(2,1),xTraj(3,1),'color',color{mod(i-1,n_color)+1},'Linewidth',1.5);
%     scatter3(xTraj(1,end),xTraj(2,end),xTraj(3,end),'b*','Linewidth',1.5);
    figure(3);  
    subplot(2,1,1); hold on;
    h_plt_u1 = plot(t_vec(1:end-1),uTraj(1,:),'color',color{mod(i-1,n_color)+1},'Linewidth',1);
    subplot(2,1,2); hold on;
    h_plt_u2 = plot(t_vec(1:end-1),uTraj(2,:),'color',color{mod(i-1,n_color)+1},'Linewidth',1);
end
figure(2);
scatter3(0,0,0,200,'r','filled','p','MarkerFaceColor','r');
view(130.2,37.8);
xlabel('$x_1$','interpreter','latex');ylabel('$x_2$','interpreter','latex');zlabel('$x_3$','interpreter','latex');
goodplot([6 6]);
legend off;
set(gca, 'box', 'off');
if print_file 
    fig_name = ['exp2_ulim_' num2str(include_input_limits) '_min_norm_', num2str(min_norm_control) '_3d_trajs']; 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end
figure(3)
subplot(2,1,1)
ylabel('$u_1$','interpreter','latex')
plot(t_vec,-1*ones(size(t_vec)),'r--',t_vec,1*ones(size(t_vec)),'r--')
grid on;
ylim([-1.1, 1.1])
goodplot([6 4])

subplot(2,1,2)
plot(t_vec,-1*ones(size(t_vec)),'r--',t_vec,1*ones(size(t_vec)),'r--')
ylim([-1.1, 1.1])
ylabel('$u_2$','interpreter','latex')
xlabel('Time (s)');
% legend([h_plt_u1 h_plt_u2],{'u_1','u_2'},'Orientation','Horizontal','interpreter','latex');
goodplot([6 4]);
% goodplot([6 2]);
% set(gca, 'box', 'off');
grid on;%ylim([-1.1 1.1])
if print_file 
    fig_name = ['exp2_ulim_' num2str(include_input_limits) '_min_norm_', num2str(min_norm_control) '_trajs']; 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end
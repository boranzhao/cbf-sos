
options = mskoptimset;

% simulate
t_vec = 0:dt:duration;
T_steps = length(t_vec)-1;
xTraj = zeros(n,T_steps+1);
uTraj = zeros(nu,T_steps);
x = x0;
xTraj(:,1) = x;

for i=1:T_steps
    t = t_vec(i);  
    if min_norm_control == 0
        u  = u_fcn(x);           % using the control law output by the SOS optimization
    else
        % min-norm control law with the searched cbfh_fcn1
        dh_dx0 = dh_dx_fcn(x);
        % phi0 + phi1*u >=0 
        phi0 = dh_dx0*f_fcn(x) + cbf_config.alpha*h_fcn(x);
        phi1 = dh_dx0*B;  
        % without constraint
%         if phi0>=0
%             u = zeros(nu,1);
%         else
%             u = -phi1'*phi0/(phi1*phi1');
%         end 
        
        % additionally add constraint            
        [u,feval,exitflag,~] = quadprog(eye(nu),0,-phi1,phi0,[],[],-1./diag(Du),1./diag(Du),[],options); % options
    end
    % record
    uTraj(:,i) = u;        
    % propagate with zero-hold for control inputs
    [d_t,d_state] = ode23(@(t,x) dynamics(t,x,u),[t_vec(i) t_vec(i+1)],x); %ode_options
    
    % update and record;
    x = d_state(end,:)';    
    xTraj(:,i+1) = x;    
end  

% plot the result
figure(1)
plot(xTraj(1,:),xTraj(2,:),'k-.','linewidth',1);
scatter(xTraj(1,1),xTraj(2,1),'ko')
legend([plt_hinit plt_h],{'$h^\textrm{init}(x)=0$','$h(x)=0$'},'interpreter','latex','Orientation','Horizontal'); %'Orientation','Horizontal'

goodplot([3.5 4]);
fig_name = ['exp1_ulim_' num2str(include_input_limits) '_linearlike_' num2str(linear_like_form) '_cbf2'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end	

figure(2);hold on;
plot(t_vec(1:end-1),uTraj);
plot(t_vec,-1*ones(size(t_vec)),'r--',t_vec,1*ones(size(t_vec)),'r--')
ylabel('$u$','interpreter','latex');
xlabel('Time (s)');
ylim([-1.1, 1.1])
goodplot([6 3])
fig_name = ['exp1_ulim_' num2str(include_input_limits) '_linearlike_' num2str(linear_like_form) '_u'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end	

if save_file 
    file_name = ['exp1_ulim_' num2str(include_input_limits) '_linearlike_' num2str(linear_like_form) '_degY_' num2str(deg_Y) '_deg_ux_refine_' num2str(deg_ux_refine) '.mat'];    
    save(file_name,'plant','cbf_config','h_init_fcn','h_init_fcn1',...
        'h_fcn','h_fcn1','u_fcn','dh_dx_fcn','X_fcn','Y_fcn','min_eig_Qhs','t_vec','xTraj','uTraj');
end
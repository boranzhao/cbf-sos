% specify the first file
load exp1_ulim_1_linearlike_1_degY_2_deg_ux_refine_3

figure;hold;
plot(t_vec(1:end-1),uTraj,'k','linewidth',1);

% specific the second file
load exp1_ulim_1_linearlike_2_degY_2_deg_ux_refine_5

plot(t_vec(1:end-1),uTraj,'k-.','linewidth',1);
plot(t_vec,-1*ones(size(t_vec)),'r--',t_vec,1*ones(size(t_vec)),'r--')

ylabel('$u$','interpreter','latex');
xlabel('Time (s)');
ylim([-1.1, 1.1])

legend('$z(x)=x$','$z(x)=[x_1\ x_2\ x_2^3]^T$','interpreter','latex')
goodplot([6,2])
fig_name = ['exp1_ulim_' num2str(include_input_limits) '_u'];
if print_file 
    savefig([fig_name '.fig']);
    print([fig_name '.pdf'], '-painters', '-dpdf', '-r150');
end	
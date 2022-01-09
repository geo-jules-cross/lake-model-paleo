%% --------------------Plotting---------------------

figure(1)
hold on
plot(t_vec,all_data.h_LB,'Color','r','LineWidth',2) %LB
plot(t_vec,all_data.h_LH,'g','LineWidth',2)
plot(t_vec,all_data.h_LF,'b','LineWidth',2)
line(t_vec,FH_spillpoint,'Color','k')
line(t_vec,HB_spillpoint,'Color','k')
legend('LB','LH','LF')
ylabel('m asl')
xlabel('time (year)')
hold off

% figure(2),plot(t_vec, all_data.inflow_new_LF,'Color','b')
% hold on
% plot(t_vec, all_data.inflow_new_LH,'Color','g')
% title('Inflow test')
% ylabel('Inflow m^3')
% xlabel('time (year)')
% hold off

% figure(3)
% [AX,H1,H2] = plotyy(t_vec,all_data.e_LF,t_vec,all_data.inflow_new_LF,'plot');
% set(get(AX(1),'Ylabel'),'String','LF lake level (m asl)') 
% set(get(AX(2),'Ylabel'),'String','Inflow (m^3)') 
% xlabel('time (year)')
% set(H1,'LineWidth',2)
% set(H2,'LineStyle','--')
% legend('LF elevation','Inflow')

% figure(4)
% [AX,H1,H2] = plotyy(t_vec,all_data.e_LH,t_vec,all_data.inflow_new_LH,'plot');
% set(get(AX(1),'Ylabel'),'String','LH lake level (m asl)') 
% set(get(AX(2),'Ylabel'),'String','Inflow (m^3)') 
% xlabel('time (year)')
% set(H1,'LineWidth',2)
% set(H2,'LineStyle','--')
% legend('LH elevation','Inflow')

% figure(5); subplot(2,1,1)
% plot(t_vec,all_data.e_LB,'Color','r','LineWidth',2) %LB
% hold on
% plot(t_vec,all_data.e_LH,'g','LineWidth',2)
% plot(t_vec,all_data.e_LF,'b','LineWidth',2)
% line(t_vec,LFspill,'Color','b')
% line(t_vec,LHspill,'Color','g')
% line(t_vec,LBspill,'Color','r')
% legend('LB','LH','LF')
% ylabel('m asl')
% grid on
% title('Merging of TV lakes')
% xlabel('Time (year)')
% hold off
% subplot(2,1,2)
% [AX,H1,H2] = plotyy(t_vec(1:end-1),all_data.rate_LF,t_vec,all_data.inflow_new_LF,'plot');
% set(get(AX(1),'Ylabel'),'String','LF inflow (de/dt)') 
% set(get(AX(2),'Ylabel'),'String','LF inflow (m^3)') 
% xlabel('Time (year)')
% set(H1,'LineWidth',2)
% set(H2,'LineStyle','--')
% grid on

plot(hyps_LF(:,2)/10e6,hyps_LF(:,1))
hold on
plot(hyps_LH(:,2)/10e6,hyps_LH(:,1))
plot(hyps_LB(:,2)/10e6,hyps_LB(:,1))
hold off
ylabel('Elevation (m asl)')
xlabel('Surface area (x10^6 m^2)')

%figure is exported as 6 x 3 inches
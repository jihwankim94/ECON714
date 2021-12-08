%figure8.m

clear all
clc

show_horz=18;
time=0:0.25:show_horz/4-0.25;

load simulation_results

figure;

subplot(2,3,1)

plot(time,NaN*unemp.case1b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(2:show_horz+1),'b'); hold on
plot(time,unemp.case1b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);
ylabel('Estimated Price Stickiness (\xi=0.75)');

axis([0 4 3.5 7]);

leg=legend('Benefits AR(1)=0.90','Benefits AR(1)=0.75');
set(leg,'Position',[0.2 0.6 0.07 0.04]);
legend boxoff;

subplot(2,3,4)
plot(time,unemp.case2a(2:show_horz+1),'b');hold on
plot(time,unemp.case2b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);
ylabel('More Flexible Prices (\xi=0.5)');
xlabel('Years');

axis([0 4 3.5 7]);

subplot(2,3,2)

plot(time,unemp.case3a(2:show_horz+1),'b'); hold on
plot(time,unemp.case3b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);

axis([0 4 3.5 7]);

subplot(2,3,5)

plot(time,unemp.case4a(2:show_horz+1),'b');hold on
plot(time,unemp.case4b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);
xlabel('Years');

axis([0 4 3.5 7]);

subplot(2,3,3)

plot(time,unemp.case5a(2:show_horz+1),'b'); hold on
plot(time,unemp.case5b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);

axis([0 4 3.5 7]);

subplot(2,3,6)

plot(time,unemp.case6a(2:show_horz+1),'b');hold on
plot(time,unemp.case6b(2:show_horz+1),'r-.'); hold on
plot(time,unemp.case1a(1)*ones(1,show_horz),'k--'); hold on

title('Unemployment Rate (%)','FontSize',12);
xlabel('Years');

axis([0 4 3.5 7]);

text(-9,12.5,'Normal Times','FontSize',12);
text(-4,12.5,'1 Year ZLB','FontSize',12);
text(1,12.5,'2 Years ZLB','FontSize',12);

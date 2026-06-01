% Step 1: Introduction of Problem Parameters and Initial Conditions

k_01 = 0.15; %0.15;
k_02 = 0.10; %0.10;
k_03 = 0.05; %0.05;

s_0  = 0.5;
e_0  = 0.2;
c_0  = 0.0;
p_0  = 0.2;

% Step 2: Initialization of Time Grid and Solution Components

T    = 10000;
h    = 10;
t    = (0:h:T)';

s    = zeros(length(t),1);
e    = zeros(length(t),1);
c    = zeros(length(t),1);
p    = zeros(length(t),1);

s(1) = s_0;
e(1) = e_0;
c(1) = c_0;
p(1) = p_0;

% Step 3: Explicit Reformulation of Implicit Eulerian Time-Stepping Method

for j = 1:1:(length(t)-1)
   h      = t(j+1) - t(j);
   D      = (1/(h*k_01)) + ((k_02+k_03)/(k_01)) + (e(1)+c(1)-s(j)-c(j)) + h*k_03*(e(1)+c(1));
   F      = -((1/(h*k_01)) + ((k_02+k_03)/(k_01)))*s(j) - (k_02/k_01)*c(j);
   s(j+1) = -F/((D/2) + sqrt((D/2)^2-F));
   c(j+1) = (c(j)+h*k_01*(e(1)+c(1))*s(j+1))/(1+h*k_01*s(j+1)+h*(k_02+k_03));
   e(j+1) = e(1)+c(1)-c(j+1);
   p(j+1) = s(1)+c(1)+p(1)-s(j+1)-c(j+1);
endfor

% Step 4: Calculation of largest eigenvalue of system matrix

lambda = -k_01*k_03*(e_0+c_0)/(k_01*(e_0+c_0)+k_01*(s_0+c_0+p_0)+(k_02+k_03));

% Step 5: Calculation of limit (c/c+s) =: B for t to infinity

A = - ((k_02+k_03-k_01*(e_0+c_0))/(2*k_02)) + sqrt(((k_02+k_03-k_01*(e_0+c_0))/(2*k_02))^2+((k_01*(e_0+c_0))/(k_02)));
B = A/(1+A);

% Step 5.1: alpha of Bersani 2012

alpha_bersani = (k_01/2)*(((k_02+k_03)/k_01)+(e_0+c_0))*(1-sqrt(1-4*k_03*(e_0+c_0)/(k_01*(((k_02+k_03)/k_01)+(e_0+c_0)))));
K_W_bersani   = (k_03-alpha_bersani)/alpha_bersani;

% Step 6: Plotting of Solution Components and Conservation Laws

figure(1)
plot(t,s,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,e,'linewidth',0.8,'marker','+','linestyle','-')
hold on
plot(t,c,'linewidth',0.8,'marker','square','linestyle','-')
hold on
plot(t,p,'linewidth',0.8,'marker','diamond','linestyle','-')
hold on
plot(t,(e_0+c_0)*ones(length(t),1),'linewidth',0.8,'linestyle','--')
hold on
plot(t,(s_0+c_0+p_0)*ones(length(t),1),'linewidth',0.8,'linestyle','-.')
title('Concentrations of enzymatic reaction','fontsize',20)
legend({'s(t)','e(t)','c(t)','p(t)','e^{*}','p^{*}'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('Concentrations','fontsize',14)
xlim([0 1000])
ylim([0 0.8])
hold off

##figure(2)
##plot(t,e+c,'linewidth',0.8,'marker','*')
##hold on
##plot(t,s+c+p,'linewidth',0.8,'marker','+')
##title('Conservation laws of enzymatic reaction','fontsize',20)
##legend({'e(t)+c(t)','s(t)+c(t)+p(t)'},'location','eastoutside','fontsize',14)
##xlabel('Time','fontsize',14)
##ylabel('Conservation laws','fontsize',14)
##xlim([0 1000])
##ylim([0 0.8])
##hold off

figure(3)
plot(t,s,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |s(t) - 0|','fontsize',20)
legend({'s(t)','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('s(t)','fontsize',14)
xlim([0 1000])
ylim([0 0.8])
hold off

figure(4)
plot(t,c,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |c(t) - 0|','fontsize',20)
legend({'c(t)','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('c(t)','fontsize',14)
xlim([0 1000])
ylim([0 0.8])
hold off

figure(5)
plot(t,abs(e-(e_0+c_0)*ones(length(e),1)),'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |e(t) - e^{*}|','fontsize',20)
legend({'|e(t) - e^{*}|','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('|e(t) - e^{*}|','fontsize',14)
xlim([0 1000])
ylim([0 0.8])
hold off

figure(6)
plot(t,abs(p-(s_0+c_0+p_0)*ones(length(e),1)),'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |p(t) - p^{*}|','fontsize',20)
legend({'|p(t) - p^{*}|','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('|p(t) - p^{*}|','fontsize',14)
xlim([0 1000])
ylim([0 0.8])
hold off

figure(7)
plot(t,s,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |s(t) - 0| (Zoom in)','fontsize',20)
legend({'s(t)','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('s(t)','fontsize',14)
xlim([2000 10000])
ylim([0 5E-08])
hold off

figure(8)
plot(t,c,'linewidth',0.8,'marker','*','linestyle','-')
hold on
plot(t,1*(s_0+e_0+c_0+p_0).*exp(-abs(lambda).*t),'linewidth',0.8,'linestyle','-')
title('Exponential decay of |c(t) - 0| (Zoom in)','fontsize',20)
legend({'c(t)','exp. decay'},'location','eastoutside','fontsize',14)
xlabel('Time','fontsize',14)
ylabel('c(t)','fontsize',14)
xlim([2000 10000])
ylim([0 1E-08])
hold off

figure(9)
plot(t,c./s,'color','blue','linewidth',0.8,'linestyle','-')
hold on
plot(t,A*ones(length(t),1),'color','black','linewidth',0.8,'linestyle','--')
hold off
legend({'c(t)/s(t)','Limit A'},'location','eastoutside','fontsize',14)
title('Convergence of c(t)/s(t) to A','fontsize',20)
xlabel('Time','fontsize',14)
ylabel('c(t)/s(t)','fontsize',14)
hold off

figure(10)
plot(t(2:1:end),(e(2:1:end).*s(2:1:end))./c(2:1:end),'color','blue','linewidth',0.8,'linestyle','-')
hold on
plot(t(2:1:end),((e_0+c_0)/A)*ones(length(t(2:1:end)),1),'color','black','linewidth',0.8,'linestyle','--')
hold off
legend({'(e(t)*s(t))/c(t)','Limit'},'location','eastoutside','fontsize',14)
title('Convergence of (e(t)*s(t))/c(t)','fontsize',20)
xlabel('Time','fontsize',14)
ylabel('(e(t)*s(t))/c(t)','fontsize',14)
hold off

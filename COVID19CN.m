clear; clc; clf;
delete(findall(0,'Type','figure'))

%% Generating the WS complex network
k_measure = [4 5 6 7 8 9 10 11 12 13 14]; m = 4; K = length(k_measure);

P = zeros(1,K);
for k = 1:K
    P(k) = m^(k_measure(k)-m)/factorial(k_measure(k)-m)*exp(-m);
end
for k = m:20
    Q(k) = m^(k-m)/factorial(k-m)*exp(-m);
end

figure
subplot(2,1,1)
plot(k_measure, P, 'color',[0.2 0.5 0.9],'LineWidth',2); hold on;
set(gca,'FontSize',14);
title('Watts-Strogatz Discrete Network');
xlabel('Connectivity measure','FontSize',14);
ylabel('Probability distribution','FontSize',14);
grid on; hold off;

subplot(2,1,2)
plot(Q, 'color',[0.9 0.5 0.2],'LineWidth',2); hold on;
set(gca,'FontSize',14);
title('Watts-Strogatz Full Plot');
xlabel('Connectivity measure','FontSize',14);
ylabel('Probability distribution','FontSize',14);
grid on; hold off;

%%  Parameters
% setting general parameters for simulations
T = 1000; dt = 0.1275; % time horizon & time step
gamma = .447; lambda = .469; % gamma and lambda in SAIR
sigma = .087; mu = .0523; % sigma and mu in SAIR
alpha = 0.115;% alpha in SAIR
% g = 0.458345798579201676; l = g+0.0191; % gamma and lambda in SAIR
% s = 0; m = 0.052857235642550742; % sigma and mu in SAIR
% a = 0.11118405987804617;% alpha in SAIR

N = 10^7;
%k_measure = [2 4 6 8 10 12]; K = length(k_measure);
k_classes = [0.0183 0.0741 0.1469 0.1957 0.1957 0.1566 0.1045 0.0597 0.0299 0.0133 0.0053];% 0.000221];

psi_1 = k_measure./20%[0.05 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7];
psi_2 = (k_measure./20)*0.4%[0.05 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7].*0.4;

% psi_1 = ones(1,6);
% psi_2 = ones(1,6);

%% DYNAMICS
S = zeros(K,T); A = zeros(K,T); I = zeros(K,T); R = zeros(K,T);
dS = zeros(K,T); dA = zeros(K,T); dI = zeros(K,T); dR = zeros(K,T);
dDR = zeros(K,T); % increment of diagnosed removed
Atot = zeros(K,T); %Atot(1,1) = A(1,1); Atot(2,1) = A(2,1);
Itot = zeros(K,T); %Itot(1,1) = I(1,1); Itot(2,1) = I(2,1);
%f = sigma.*ones(T,2);
for k=1:K
    A(k,1) = 40/(N*k_classes(k));
    I(k,1) = 5/(N*k_classes(k));
    S(k,1) = 1 - A(k,1) - I(k,1);
    Atot(k,1) = A(k,1);
    Itot(k,1) = I(k,1);
end



% initialising theta_1 and theta_2
theta_1 = zeros(1,T); theta_2 = zeros(1,T);
for k = 1:K
        Pk = k_measure(k)*k_classes(k);
    
        theta_1(1) = theta_1(1) + Pk*A(k,1);
        theta_2(1) = theta_2(1) + Pk*I(k,1);
end
theta_1(1) = theta_1(1)/(2*m);
theta_2(1) = theta_2(1)/(2*m);


for t = 2:T
    
    for k = 1:K
        
        dS(k,t) = (-S(k,t-1)*(psi_1(k)*gamma*theta_1(t-1) + psi_2(k)*lambda*theta_2(t-1)))*dt;
        dA(k,t) = (S(k,t-1)*(psi_1(k)*gamma*theta_1(t-1) + psi_2(k)*lambda*theta_2(t-1))-A(k,t-1)*(sigma+alpha))*dt;
        dI(k,t) = (alpha*A(k,t-1)-mu*I(k,t-1))*dt;

        %dDR(t) = (mu*I(k,t-1))*dt;

        S(k,t) = S(k,t-1) + dS(k,t);
        A(k,t) = A(k,t-1) + dA(k,t);
        I(k,t) = I(k,t-1) + dI(k,t);
        
        Atot(k,t) = Atot(k,t-1) - dS(k,t);
        Itot(k,t) = Itot(k,t-1) + alpha*A(k,t-1)*dt;

        %dR(t) = S(t) + A(t) + I(t);
        R(k,t) = 1 - S(k,t) - A(k,t) - I(k,t);
        
        % calculating theta_1 and theta_2
        
        Pk = k_measure(k)*k_classes(k);
    
        theta_1(t) = theta_1(t) + Pk*A(k,t);
        theta_2(t) = theta_2(t) + Pk*I(k,t);

    end
    theta_1(t) = theta_1(t)/(2*m);
    theta_2(t) = theta_2(t)/(2*m);
    
    if t > 10
        psi_1 = psi_1-(psi_1.*0.005); psi_2 = psi_2-(psi_2.*0.00135);
    end
end


% syms x1 x2 x3 g l s m
% [Sx1,Sx2,Sx3] = solve(-x1*(g*x2+l*x3) == 0, g*x1*x2 - s*x2 == 0, l*x1*x3 - m*x3 ==0)

syms t1 t2 k1 k2 g l s m a
X = [-k1*gamma*t1-k2*lambda*t2 0 0 0;
    -k1*gamma*t1-k2*lambda*t2 -alpha-sigma 0 0;
    0 alpha -mu 0;
    0 sigma mu 0];
det(X)



% ====================================================================
% =======================GRADE/POPULATION=============================
% plot(S(1,:),k(1,:),'--','color','r'); hold on;
% plot(S(8,:),k(1,:),'--','color','b'); hold on;
% plot(S(13,:),k(1,:),'-','color','b'); hold on;
% plot(S(18,:),k(1,:),'--','color','g'); hold on;
% plot(S(20,:),k(1,:),'-','color','g')
% title('Evolutionary Gamk(2,kk)e Dynamics of Honeybees Swarm')
% xlabel('Population %')
% ylabel('Grade of the nodes')
% legend('S1','S8','S13','S18','S20')

%% =======================POPULATION/TIME==============================
% cAI = zeros(1,T);
% for t=1:T
%     cAI(t) = sum(dA(1:t))+sum(dI(1:t));
% end

figure
t = datetime(2020,1,22) + caldays(0:T-1);
colormap jet(11)
jetcustom = jet(11);
subplot(2,1,1)
for k=1:K
    plot(1:T,Atot(k,1:T).*(N*k_classes(k)),'-','Color',jetcustom(k,:),'LineWidth',2); hold on; %A(1:T).*(58.5*10^6)+
end
ylim([0 1000]);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = 0:100:T;
ax.YAxis.TickValues = 0:200:1000;
title('Watts-Strogatz Model: Asymptomatic Infected');
xlabel('Time (days)','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
colormap jet(11)
cb = colorbar; caxis([3 14]);
cb.Label.String = 'Node degree';
grid on; hold off;

subplot(2,1,2)
for k=1:K
    plot(1:T,Itot(k,1:T).*(N*k_classes(k)),'-','Color',jetcustom(k,:),'LineWidth',2); hold on; %A(1:T).*(58.5*10^6)+
end
ylim([0 1000]);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = 0:100:T;
ax.YAxis.TickValues = 0:200:1000;
title('Watts-Strogatz Model: Symptomatic Infected');
xlabel('Time (days)','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
colormap jet(11)
cb = colorbar; caxis([3 14]);
cb.Label.String = 'Node degree';
grid on; hold off;

figure
t = datetime(2020,1,22) + caldays(0:T-1);
colormap jet(11)
jetcustom = jet(11);
for k=1:K
    plot(1:T,R(k,1:T).*(N*k_classes(k)),'-','Color',jetcustom(k,:),'LineWidth',2); hold on;
end
ylim([0 1000]);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = 0:100:T;
ax.YAxis.TickValues = 0:200:1000;
title('Watts-Strogatz Model: Removed');
xlabel('Time (days)','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
colormap jet(11)
cb = colorbar; caxis([3 14]);
cb.Label.String = 'Node degree';
grid on; hold off;
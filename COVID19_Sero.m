clear; clc; clf;
delete(findall(0,'Type','figure'))

%% DATA
%A = csvread('dpc_nazionale.csv');

%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/leonardo/Google Drive/PhD/MATLAB/Research/dpc_nazionale.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2020/08/11 12:36:09

% Initialize variables.
filename = '/Users/leonardostella/Google Drive/PhD/MATLAB/Research/COVID-19 SIAM/dpc_nazionale.csv';
delimiter = ',';
startRow = 2;

% Format for each line of text:
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: text (%s)
%   column13: text (%s)
%	column14: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%f%f%f%f%f%f%f%f%s%s%f%*s%*s%*s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% Create output variable
dpcnazionale = table(dataArray{1:end-1}, 'VariableNames', {'ricoverati_con_sintomi','terapia_intensiva','totale_ospedalizzati','isolamento_domiciliare','totale_positivi','variazione_totale_positivi','nuovi_positivi','dimessi_guariti','deceduti','casi_da_sospetto_diagnostico','casi_da_screening','totale_casi'});

% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

asymptomatic_isolated = table2array(dpcnazionale(:,4));
symptomatic_hospitalised = table2array(dpcnazionale(:,3));

cumulative_infected = table2array(dpcnazionale(:,12));
cumulative_recovered = table2array(dpcnazionale(:,8));
cumulative_deaths = table2array(dpcnazionale(:,9));


%%  Parameters
% setting general parameters for simulations
T = length(cumulative_infected); dt = 1; % time horizon & time step
g = 0.46952; l = 0.48521; % gamma and lambda in SAIR
s = 0.065501; m = 0.15004; % sigma and mu in SAIR
a = 0.050017;% alpha in SAIR

mult = 10.5*10^3;
pop = [60*10^6 47*10^6];
k1 = [0.99209 .92]; k2 = [0.65056 .32];
k1_init = k1; k2_init = k2;

R0 = zeros(2,T);
R0(:,1) = (k1*g*m+k2*l*a)/((s+a)*m); % about 2.38


%% DYNAMICS
S = zeros(2,T); A = zeros(2,T); I = zeros(2,T); R = zeros(2,T);
A(1,1) = 94/pop(1); A(2,1) = 150/pop(2);
I(1,1) = 127/pop(1); I(2,1) = 2/pop(2);
S(1,1) = 1 - A(1,1) - I(1,1); S(2,1) = 1 - A(2,1) - I(2,1);

% model equations on the grade k of the node
dS = zeros(2,T); dA = zeros(2,T); dI = zeros(2,T); dR = zeros(2,T);
dDR = zeros(2,T); % increment of diagnosed removed
Atot = zeros(2,T); Atot(1,1) = A(1,1); Atot(2,1) = A(2,1);
Itot = zeros(2,T); Itot(1,1) = I(1,1); Itot(2,1) = I(2,1);
%f = sigma.*ones(T,2);
for t = 2:T
    %f(t,i) = (sigma(i)+(sigma(i)/2.5)*sin(t*sigma(i)));
    if t == 10 % after day 30, lower parameters due to the italian government's recommendations
%         %g = .441; l = .443;
%         k1 = [.65 1]; k2 = [.25 1];
%         R0 = (k1*g*m+k2*l*a)/((s+a)*m)
%     elseif t == 14 % after day 30, lower parameters due to the italian government's recommendations
%         %g = .441; l = .443;
%         k1 = [.75 1]; k2 = [.24 1];
%         R0 = (k1*g*m+k2*l*a)/((s+a)*m)
    elseif t >= 10 && t <= 12
        k1 = k1 - k1*0.09150; 
        k2 = k2 - k2*0.30000;
       
        if t == 14
            %R0(:,t) = (k1*g*m+k2*l*a)/((s+a)*m);
        end
        %R0 = (g*k1+l*k2)/(s+m)
    elseif t > 12 && t < 50
        k1 = k1 - k1*0.03850; 
        k2 = k2 - k2*0.05000;
        
    elseif t == 160
        %k1
        %k2
        
    elseif t > 160 && t < 185
        k1 = k1 + k1*0.01500; 
        k2 = k2 + k2*0.03000;
        %R0 = (g*k1+l*k2)/(s+m)    
    end
    
    for i = 1:2
        
        dS(i,t) = (-S(i,t-1)*(k1(i)*g*A(i,t-1) + k2(i)*l*I(i,t-1)))*dt;
        dA(i,t) = (S(i,t-1)*(k1(i)*g*A(i,t-1) + k2(i)*l*I(i,t-1))-A(i,t-1)*(s+a))*dt;
        dI(i,t) = (a*A(i,t-1)-m*I(i,t-1))*dt;

        dDR(t) = (m*I(i,t-1))*dt;
        
%        real_recovered(t) = real_recovered(t-1) + m*I(i,t-1);

        S(i,t) = S(i,t-1) + dS(i,t);
        A(i,t) = A(i,t-1) + dA(i,t);
        I(i,t) = I(i,t-1) + dI(i,t);
        
        Atot(i,t) = Atot(i,t-1) - dS(i,t);
        Itot(i,t) = Itot(i,t-1) + a*A(i,t-1)*dt;

        %dR(t) = S(t) + A(t) + I(t);
        R(i,t) = 1 - S(i,t) - A(i,t) - I(i,t);
    end
    
    R0(:,t) = (k1*g*m+k2*l*a)/((s+a)*m);
end

% syms x1 x2 x3 g l s m
% [Sx1,Sx2,Sx3] = solve(-x1*(g*x2+l*x3) == 0, g*x1*x2 - s*x2 == 0, l*x1*x3 - m*x3 ==0)

% syms x1 k1 k2 g l s m a
% X = [0 -k1*g*x1 -k2*l*x1 0;
%     0 k1*g*x1-a-s k2*l*x1 0;
%     0 a -m 0;
%     0 s m 0];
% det(X)




%% Calculation of Itilde
tbar = 50; % The instant where parameters do not change anymore
S0 = S(1,tbar); A0 = A(1,tbar); I0 = I(1,tbar); % Initial conditions at tbar
Sbar = S0.*0.9999; R0t = R0(1,tbar);

SS = Sbar; % Calculating parts of equation (2.13)
for i=1:1000 % Fixed point for (2.13)
    if tbar == 1
        logS = R0t*(S0-SS) + k2_init(1)*l/m*I0 + R0t*A0;
    else
        logS = R0t*(S0-SS) + k2(1)*l/m*I0 + R0t*A0;
    end
    SS = exp(-logS)*S0;
end

h = a+s-m;

if tbar == 1
    Itilde = 1/(2*k2_init(1)*l*SS)*(h-k1_init(1)*g*SS + ((h-k1_init(1)*g*SS)^2 + 4*a*k2_init(1)*l*SS)^(1/2));
    I(1,T)/A(1,T)
else
    Itilde = 1/(2*k2(1)*l*SS)*(h-k1(1)*g*SS + ((h-k1(1)*g*SS)^2 + 4*a*k2(1)*l*SS)^(1/2));
    I(1,T)/A(1,T)
end

Itilde2 = 1/(2*k2(1)*l*S(1,T))*(h-k1(1)*g*S(1,T) + ((h-k1(1)*g*S(1,T))^2 + 4*a*k2(1)*l*S(1,T))^(1/2));

%R(1,T).*pop(1) - cumulative_deaths(T)

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
t = datetime(2020,2,24) + caldays(0:T-1);
subplot(2,1,1)
plot(t(1:T),I(1,1:T).*pop(1),'-','color',[0.2 0.5 0.9],'LineWidth',3); hold on; %A(1:T).*(58.5*10^6)+
stem(t(1:T),symptomatic_hospitalised(1:T),'color',[0.36 0.56 0.86]);
%stem(t(1:T),cumulative_infected(1:T),'color',[0.99 0.33 0]);
xl = xline(datetime(2020,3,9),'-.','Lockdown','color','black','LineWidth',2);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xl.FontSize = 14;
%stem(t,confirmed_spain,'color','g');
%plot(1:T,confirmed_hubei,'r:^');
xlim([datetime(2020,2,24) datetime(2020,8,28)]);
ylim([0 2*10^5]);
xtickangle(45);
xtickformat('MMM-dd');
%set(gca,'XTick', 1:size(t),'XTickLabel',t,'FontSize',14);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = t(1:14:T);
ax.YAxis.TickValues = 0:0.4*10^5:2*10^5;
title('Confirmed Infected: Model vs. Data');
xlabel('Date','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
h = legend('Symptomatic','Data\_Hospitalized');
grid on; hold off;

subplot(2,1,2)
plot(t(1:T),A(1,1:T).*pop(1),'-','color',[0.5 0.9 0.2],'LineWidth',3); hold on; %A(1:T).*(58.5*10^6)+
stem(t(1:T),asymptomatic_isolated(1:T),'color',[0.56 0.86 0.36]);
xl = xline(datetime(2020,3,9),'-.','Lockdown','color','black','LineWidth',2);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xl.FontSize = 14;
%stem(t,confirmed_spain,'color','g');
%plot(1:T,confirmed_hubei,'r:^');
xlim([datetime(2020,2,24) datetime(2020,8,28)]);
ylim([0 2*10^5]);
xtickangle(45);
xtickformat('MMM-dd');
%set(gca,'XTick', 1:size(t),'XTickLabel',t,'FontSize',14);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = t(1:14:T);
ax.YAxis.TickValues = 0:0.4*10^5:2*10^5;
%title('Confirmed Infected: Model vs. Data');
xlabel('Date','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
h = legend('Aymptomatic','Data\_Isolated');
grid on; hold off;

figure
t = datetime(2020,2,24) + caldays(0:T-1);
plot(t(1:T),R(1,1:T).*pop(1),'-','color',[0.99 0.3 0.3],'LineWidth',3); hold on;
stem(t(1:T),cumulative_recovered(1:T)+cumulative_deaths(1:T),'color',[0.99 0.33 0]);
xl = xline(datetime(2020,3,9),'-.','Lockdown','color','black','LineWidth',2);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xl.FontSize = 14;
%plot(1:T,confirmed_hubei,'r:^');
xlim([datetime(2020,2,24) datetime(2020,8,28)]);
ylim([0 2*10^6]);
xtickangle(45);
xtickformat('MMM-dd');
%set(gca,'XTick', 1:size(t),'XTickLabel',t,'FontSize',14);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = t(1:14:T);
title('Recovered and Deaths: Model vs. Data');
xlabel('Date','FontSize',16);
ylabel('Number of removed','FontSize',16);
set(gca,'DefaultTextFontSize',14);
h = legend('Removed','Data');
grid on; hold off;


%% SEROPREVALENCE
figure

% Fitting the seroprevalence data with exponential
x=(1:T)';
y=(cumulative_infected.*6);
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(1,T,T);

% Creating the actual plot
t = datetime(2020,2,24) + caldays(0:T-1);
plot(t(1:T),Itot(1,1:T).*pop(1),'-','color',[0.2 0.5 0.9],'LineWidth',3); hold on;
plot(t(1:T),Atot(1,1:T).*pop(1),'-','color',[0.56 0.86 0.36],'LineWidth',3);
plot(t(1:T),Atot(1,1:T).*pop(1)+Itot(1,1:T).*pop(1),':','color',[0.69 0.3 0.69],'LineWidth',3);
stem(t(1:T),cumulative_infected(1:T),'color',[0.69 0.3 0.69]);
plot(t(1:T),f0(xx),'-.','color',[0.99 0.3 0.3],'LineWidth',2); hold on;
xl = xline(datetime(2020,3,9),'-.','Lockdown','color','black','LineWidth',2);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xl.FontSize = 14;
%plot(1:T,confirmed_hubei,'r:^');
xlim([datetime(2020,2,24) datetime(2020,8,28)]);
ylim([0 2*10^6]);
xtickangle(45);
xtickformat('MMM-dd');
%set(gca,'XTick', 1:size(t),'XTickLabel',t,'FontSize',14);
set(gca,'FontSize',14);
ax = gca;
ax.XAxis.TickValues = t(1:14:T);
title('Analysis: Seroprevalence Study vs Data');
xlabel('Date','FontSize',16);
ylabel('Number of cases','FontSize',16);
set(gca,'DefaultTextFontSize',14);
h = legend('Model: Symptomatic','Model: Asymptomatic','Model: Total Infected','Data\_Cumulative\_Infected','Data\_Seroprevalence\_Study','Location','northwest');
grid on; hold off;
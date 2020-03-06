set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% loading data %%

pot="double";
beta=120;
number=1;

file=sprintf('%s_%d',pot,beta)
output=sprintf('../simulations/%s_%d/%s',file,number,file)
output="simulations/output"

data=load(output+"_pot.out");
xx=data(:,1);
VV=data(:,2);
figure 
plot(xx,VV)
n_part=3;

data=load(output+"_pos.out");
n_MCS=size(data,1)
n_slices=size(data,2)/n_part
A=zeros(n_MCS,n_slices,n_part);

for k=1:n_part
	A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
end

% moy=mean(A(:,:,1),2)
% for k=1:n_slices
%     A(:,k,:)=A(:,k,:)-moy;
% end

B=zeros(n_MCS*n_slices,n_part);
for k=1:n_part
    for l=1:n_MCS
        B((1+(l-1)*n_slices):l*n_slices,k)=A(l,:,k);
    end
end
sa=mean(B(:,2))
B=B-sa;

A2=zeros(n_MCS,n_slices,n_part);
moy=A(:,:,1);
for k=1:n_part
    A2(:,:,k)=A(:,:,k)-moy;
end


B2=zeros(n_MCS*n_slices,n_part);
for k=1:n_part
    for l=1:n_MCS
        B2((1+(l-1)*n_slices):l*n_slices,k)=A2(l,:,k);
    end
end

%%
% X2=mean(A.^2,2);
% figure
% plot(1:n_MCS,X2)

% G=zeros(n_MCS,n_part);
% for k=1:n_part
%     for n=1:n_MCS
%         for i=

%% imaginary path %%

figure
k=3; % number of the particle that we want to analyse
hold on;
for k=1:n_part
    t(k)=plot([A(1,:,k) A(1,1,k)],[0:1:size(A,2)]);
end
% xlim([-10 10])
xlabel('x')
ylabel('\tau / \delta \tau')
grid on

%{
closw=waitforbuttonpress;


for i=2:n_MCS
		closw=waitforbuttonpress;
		if ~ishandle(t)
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'XData',[A(i,:,k) A(i,1,k)])
        end
end

%}

% fig3=figure('Position',[50,50,650,450]);
% hold on;
% for k=[3001 5001 7000]
%     t(k)=plot([A(k,:,1) A(k,1,1)],[0:1:size(A,2)],'LineWidth',2);
% end
% hold off
% xlim([-3 3])
% xlabel('$\tilde{x}$')
% ylabel('$\tau / \delta \tau$')
% yticks(0:16)
% yticklabels({'0' '' '2' '' '4' '' '6' '' '8' '' '10' '' '12' '' '14' '' '16'})
% xticks(-3:3:3)
% grid on, box on
% print(fig3,'threePaths', '-depsc');



%}

n=1:n_slices/2;
G=twoP_corr(A,1,n,100);

%%
figure
semilogy(n,G,'+')
xlabel('$\Delta \tau$');
ylabel('G$(\Delta \tau)$');

%% 
meff=(0.5*log(G(1:end-2)./G(3:end)));
figure 
plot(n(2:end-1),meff)

xlabel('$\Delta \tau$');
ylabel('$m_{eff}(\Delta \tau)$');

%%
Oi=mean(A(:,:,1),2);
tMC=0:n_MCS-1;
ACT=zeros(size(tMC));
for i=1:n_MCS-1
    ACT(i)=1/(n_MCS-tMC(i)-1)*sum((Oi(1:n_MCS-tMC(i))-mean(Oi(1:n_MCS-tMC(i)))).*(Oi(tMC(i)+1:n_MCS)-mean(Oi(tMC(i)+1:n_MCS))));
    if ACT(i) <= 0
       iFinal=i-1;
       break;
    end
        
end
ACT0=ACT(1);
ACT_norm=ACT(2:iFinal)/ACT0;
tMC_norm=tMC(2:iFinal);
figure 
loglog(tMC_norm,ACT_norm)
xlabel('$t_{MC}$')
ylabel('$A_x$')

tau0=0.5+sum(ACT_norm)
Neff=n_MCS/(2*tau0)
%% histogram %%%{

xMin=-10;
xMax=+10;
Ninter=10000;
X=linspace(xMin,xMax,Ninter+1);
dx=(xMax-xMin)/Ninter;

m=1;
w=1;

% V=0.1*X.^4-5*X.^2;
% V=0.5*X.^2;
% x0=3;
% V0=0.5;
% V=V0*((X./x0).^2-1).^2;

T=1;
kB=1;%1.380649*10^(-23);
hbar = 1.0545718;
beta=120%1.0/(T);

V=0.5*m*w^2*X.^2;

figure
plot(X,V)

z=partitionFct(X,V,beta,dx);

p=exp(-V*beta)/z;
%%
figure
% plot(X,p)

% pTot=0.5*dx*(p(1)+p(end));
% for i=2:length(p)-1
%     pTot=pTot+dx*p(i);
% end
% disp(pTot)

mean_x=zeros(1,3);
for k=1:3
    mean_x(k)=mean(B(:,k));
end

hold on;
for k=1:n_part
    h(k)=histogram(B(:,k)-mean_x(1),51,'Normalization','pdf');
    plot([mean_x(k) mean_x(k)]-mean_x(1),[0 25])
end

figure
% plot(X,p)

% pTot=0.5*dx*(p(1)+p(end));
% for i=2:length(p)-1
%     pTot=pTot+dx*p(i);
% end
% disp(pTot)

hold on;
for k=2:n_part
    h(k)=histogram(B2(:,k),51,'Normalization','pdf');
end

%%
figure
hold on;
for k=1:n_part
    t(k)=histogram(A(1,:,k),91,'Normalization','probability');
end
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:n_MCS
	pause(.01)
		if ~ishandle(t(1))
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'Data',A(i,:,k))
        end
        
end
%}

function Z=partitionFct(X,V,beta,dx)
if ~isequal(size(X),size(V))
    disp("fff")
else
    Z=0;
    for Ei=V
        Z=Z+exp(-Ei*beta)*dx;
    end
end
end

function G=twoP_corr(A,particle,n,n_therm)
N_MC=size(A,1);
N_slice=size(A,2);
G=zeros(size(n));
for t=(n_therm+1):N_MC
   for i=1:N_slice
       G=G+A(t,i,particle)*A(t,mod(i+n-1,N_slice)+1,particle);
   end
end
G=G/(N_MC-n_therm)/N_slice-mean(mean(A(:,:,particle),1),2)^2;
end


%{



%% loading data %%

pot="double";
beta=120;
number=1;

file=sprintf('%s_%d',pot,beta)
output=sprintf('../simulations/%s_%d/%s',file,number,file)
output="simulations/output"

data=load(output+"_pot.out");
xx=data(:,1);
VV=data(:,2);
figure
plot(xx,VV)
n_part=3;

data=load(output+"_pos.out");
n_MCS=size(data,1)
n_slices=size(data,2)/n_part
n_therm=1
n_MCS=n_MCS-n_therm
A=zeros(n_MCS,n_slices,n_part);
B=zeros(n_MCS*n_slices,n_part);
for k=1:n_part
	A(:,:,k)=data(n_therm+1:end,(1+(k-1)*n_slices):k*n_slices);
    for l=1:n_MCS
        B((1+(l-1)*n_slices):l*n_slices,k)=A(l,:,k);
    end
end

%%
% X2=mean(A.^2,2);
% figure
% plot(1:n_MCS,X2)

% G=zeros(n_MCS,n_part);
% for k=1:n_part
%     for n=1:n_MCS
%         for i=

%% imaginary path %%
%{
figure
%k=3; % number of the particle that we want to analyse
hold on;
for k=1:n_part
    t(k)=plot([A(1,:,k) A(1,1,k)],[0:1:size(A,2)]);
end
xlim([-10 10])
xlabel('x')
ylabel('\tau / \delta \tau')
grid on
closw=waitforbuttonpress;


for i=2:n_MCS
		pause(.01)
		if ~ishandle(t)
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'XData',[A(i,:,k) A(i,1,k)])
        end:
end

%}

% fig3=figure('Position',[50,50,650,450]);
% hold on;
% for k=[3001 5001 7000]
%     t(k)=plot([A(k,:,1) A(k,1,1)],[0:1:size(A,2)],'LineWidth',2);
% end
% hold off
% xlim([-3 3])
% xlabel('$\tilde{x}$')
% ylabel('$\tau / \delta \tau$')
% yticks(0:16)
% yticklabels({'0' '' '2' '' '4' '' '6' '' '8' '' '10' '' '12' '' '14' '' '16'})
% xticks(-3:3:3)
% grid on, box on
% print(fig3,'threePaths', '-depsc');





n=1:n_slices/2;
G=twoP_corr(A,1,n,100);

%%
figure
semilogy(n,G,'+')
xlabel('$\Delta \tau$');
ylabel('G$(\Delta \tau)$');

%%
meff=(0.5*log(G(1:end-2)./G(3:end)));
figure
plot(n(2:end-1),meff)

xlabel('$\Delta \tau$');
ylabel('$m_{eff}(\Delta \tau)$');

%%

Oi=mean(A(:,:,1),2);
tMC=0:n_MCS-1;
ACT=zeros(size(tMC));
for i=1:n_MCS-1
    ACT(i)=1/(n_MCS-tMC(i)-1)*sum((Oi(1:n_MCS-tMC(i))-mean(Oi(1:n_MCS-tMC(i)))).*(Oi(tMC(i)+1:n_MCS)-mean(Oi(tMC(i)+1:n_MCS))));
    if ACT(i) <= 0
       iFinal=i-1;
       break;
    end

end
ACT0=ACT(1);
ACT_norm=ACT(2:iFinal)/ACT0;
tMC_norm=tMC(2:iFinal);
figure
loglog(tMC_norm,ACT_norm)
xlabel('$t_{MC}$')
ylabel('$A_x$')

tau0=0.5+sum(ACT_norm)
Neff=n_MCS/(2*tau0)
%% histogram %%%{

xMin=-10;
xMax=+10;
Ninter=10000;
X=linspace(xMin,xMax,Ninter+1);
dx=(xMax-xMin)/Ninter;

m=1;
w=1;

% V=0.1*X.^4-5*X.^2;
% V=0.5*X.^2;
% x0=3;
% V0=0.5;
% V=V0*((X./x0).^2-1).^2;

T=1;
kB=1;%1.380649*10^(-23);
hbar = 1.0545718;
beta=120%1.0/(T);

V=0.5*m*w^2*X.^2;

figure
plot(X,V)

z=partitionFct(X,V,beta,dx);

p=exp(-V*beta)/z;
%%
figure
% plot(X,p)

% pTot=0.5*dx*(p(1)+p(end));
% for i=2:length(p)-1
%     pTot=pTot+dx*p(i);
% end
% disp(pTot)

hold on;
for k=1:n_part
    h(k)=histogram(B(:,k),151,'Normalization','pdf');
end

%%
figure
hold on;
for k=1:n_part
    t(k)=histogram(A(1,:,k),91,'Normalization','probability');
end
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:n_MCS
	pause(.01)
		if ~ishandle(t(1))
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'Data',A(i,:,k))
        end

end


function Z=partitionFct(X,V,beta,dx)
if ~isequal(size(X),size(V))
    disp("fff")
else
    Z=0;
    for Ei=V
        Z=Z+exp(-Ei*beta)*dx;
    end
end
end

function G=twoP_corr(A,particle,n,n_therm)
N_MC=size(A,1);
N_slice=size(A,2);
G=zeros(size(n));
for t=(n_therm+1):N_MC
   for i=1:N_slice
       G=G+A(t,i,particle)*A(t,mod(i+n-1,N_slice)+1,particle);
   end
end
G=G/(N_MC-n_therm)/N_slice-mean(mean(A(:,:,particle),1),2)^2;
end

%}

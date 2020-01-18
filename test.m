set(groot,'defaultAxesFontSize',13)
set(groot,'defaultAxesLabelFontSizeMultiplier',19/13)
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% loading data %%
n_part=2;
data=load("output2_pos.out");
n_MCS=size(data,1);
n_slices=size(data,2)/n_part;
A=zeros(n_MCS,n_slices,n_part);
B=zeros(n_MCS*n_slices,n_part);
for k=1:n_part
	A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
    for l=1:n_MCS
        B((1+(l-1)*n_slices):l*n_slices,k)=A(l,:,k);
    end
end

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
        end
end

%}

fig3=figure('Position',[50,50,650,450]);
hold on;
for k=[3001 5001 7000]
    t(k)=plot([A(k,:,1) A(k,1,1)],[0:1:size(A,2)],'LineWidth',2);
end
hold off
xlim([-3 3])
xlabel('$\tilde{x}$')
ylabel('$\tau / \delta \tau$')
yticks(0:16)
yticklabels({'0' '' '2' '' '4' '' '6' '' '8' '' '10' '' '12' '' '14' '' '16'})
xticks(-3:3:3)
grid on, box on
print(fig3,'threePaths', '-depsc');



%}

%% histogram %%%{

xMin=-10;
xMax=+10;
Ninter=100000;
X=linspace(xMin,xMax,Ninter+1);
dx=(xMax-xMin)/Ninter;

V=0.1*X.^4-5*X.^2;
V=0.5*X.^2;
x0=3;
V0=0.5;
V=V0*((X./x0).^2-1).^2;

T=0.1;
kB=1;%1.380649*10^(-23);
beta=1.0/(T);

figure
plot(X,V)

z=partitionFct(X,V,beta,dx);

p=exp(-V*beta)/z;

figure
plot(X,p)

pTot=0.5*dx*(p(1)+p(end));
for i=2:length(p)-1
    pTot=pTot+dx*p(i);
end
disp(pTot)

hold on;
for k=1:n_part
    h(k)=histogram(B(:,k),101,'Normalization','pdf');
end

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


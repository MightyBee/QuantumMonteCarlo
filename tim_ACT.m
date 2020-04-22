set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% loading data %%

pot="harm";
beta=120;
number=1;
n_part=1;


file=sprintf('%s_%d',pot,beta)
output=sprintf('simulations/%s_%d/%s',file,number,file)



data=load(output+"_pos.out");
n_MCS=size(data,1)
n_slices=size(data,2)/n_part
A=zeros(n_MCS,n_slices,n_part);

for k=1:n_part
	A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
end


B=zeros(n_MCS*n_slices,n_part);
for k=1:n_part
    for l=1:n_MCS
        B((1+(l-1)*n_slices):l*n_slices,k)=A(l,:,k);
    end
end


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
tau0_int=0.5+sum(ACT_norm)
Neff=n_MCS/(2*tau0)
 %%
l=[1500, 20, 11];
l=l(number);
X=tMC_norm(1:l)';
Y=log(ACT_norm(1:l))';
[yFit,a]=fit_lin(X,Y,true);
figure
plot(X,Y,X,yFit)
yFit=exp(yFit);
tau0_exp=-1/a

%%
f=figure('Position',  [200, 200, 600, 400]);
semilogy(tMC_norm,ACT_norm,'-')
hold on
semilogy(X,yFit,'-')
xlabel('$t_{\rm MC}$')
ylabel('$A_x/\sigma_x^2$')
lgd=legend('Numerical data', 'Exponential fit');
set(lgd,'fontsize',14,'Location','southwest');
print(f,sprintf('FIG/act_%s_%d',file,number), '-depsc');
saveas(f,sprintf('FIG/act_%s_%d.fig',file,number));



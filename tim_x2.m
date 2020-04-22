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

for number=1:3

file=sprintf('%s_%d',pot,beta);
output=sprintf('simulations/%s_%d/%s',file,number,file);

data=load(output+"_pos.out");
n_MCS=size(data,1)
n_slices=size(data,2)/n_part
A=zeros(n_MCS,n_slices,n_part);


for k=1:n_part
	A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
end



%% analysis %%
X2=mean(A.^2,2);

%% figure %%
f=figure('Position',  [200, 200, 600, 350]);
plot((0:20:n_MCS-1)/20,X2(1:20:n_MCS),'.')
xlabel('Metropolis Sweeps')
ylabel('$\langle x^2 \rangle$')
xlim([-20 2000+20])
ylim([0 8])
rm_space(gca)
print(f,sprintf('FIG/x2_%s_%d',file,number), '-depsc');
% saveas(f,sprintf('FIG/x2_%s_%d.fig',file,number));

end
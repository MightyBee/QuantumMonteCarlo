set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultFigurePosition',  [200, 200, 450, 300]);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '';
code = 'main3';
dossier='simulations/';


xc2=2.7;

name="xc2";
param = xc2; %


%% Simulations %%
%%%%%%%%%%%%%%%%%


output = sprintf('%s%s=%.15g.out',dossier,name,param) ;
% Execution du programme en lui envoyant la valeur a scanner en argument
cmd = sprintf('./%s%s config/Hbonds.in %s=%.15g output=%s', repertoire, code, name, param,output);
disp(cmd)
% system(cmd);


%% Analyse %%
%%%%%%%%%%%%%

data=load(output+"_e0.out");
Heth=data(1,1);
errHeth=data(1,2);
data=load(output+"_pos.out");
n_part=3;
n_MCS=size(data,1)
n_slices=size(data,2)/n_part;
A=zeros(n_MCS,n_slices,n_part);

for k=1:n_part
    A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
end
X2=mean(A(:,:,2).^2,2);
n_therm=101
A_therm=A(n_therm+1:end,:,:);
n_MCS_therm=n_MCS-n_therm

B=zeros(n_MCS_therm*n_slices,n_part);
for k=1:n_part
    for l=1:n_MCS_therm
        B((1+(l-1)*n_slices):l*n_slices,k)=A_therm(l,:,k);
    end
end
sa=mean(B(:,2));
% B=B-sa;

mean_x=zeros(1,3);
std_x=zeros(1,3);
for k=1:3
    mean_x(k)=mean(B(:,k));
    std_x(k)=std(B(:,k));
end
Rmoy=mean_x(3)-mean_x(1)
% B=B-Rmoy/2;
%     Rmoy=mean(B(:,3)-B(:,1))
%     errx=sqrt(std_x(1)^2+std_x(3)^2);
errx=std(std_x(1)+std_x(3))/sqrt(n_MCS_therm*n_slices);

f=1;
fig(f)=figure('Position',[200,400,600,400]);
box on
hold on;
for k=1:n_part
    h(k)=histogram(B(:,k),-3:0.018:3,'Normalization','pdf');
    %         plot([mean_x(k) mean_x(k)]-mean_x(1),[0 25])
end
xlim([-1.5 1.5])
xlabel('$x \; {\rm[\AA]}$')
ylabel('Probability density')
rm_space(gca);
print(fig(f),sprintf('FIG/OH3part_psi2_xc2=%d',xc2*100), '-depsc');

f=f+1;
fig(f)=figure('Position',[200,400,600,400]);
plot(1:n_MCS,X2)
xlabel('$x \; {\rm[\AA]}$')
ylabel('Probability density')
rm_space(gca);
print(fig(f),'FIG/OH3part_therm', '-depsc');


function y = morse(x)
    global R;
    
    D = 83.402;
    a = 2.2;
    r0 = 0.96;
    delta1 = 0.4 * D;
    b = 2.2;
    R1 = 2 * r0 + 1/a;
    DELTA = delta1 * exp(-b*(R-R1));

    y = (0.5*((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))+(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))) - sqrt(((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))-(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0))))).^2+4*DELTA^2)));
end



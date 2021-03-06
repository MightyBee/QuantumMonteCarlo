set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '';
code = 'main2';
dossier='simulations/transition_scan/';

nsimul = 10; % number of simulations


beta=linspace(1, 15, nsimul);

% d_tau=0.1;
% N=round(beta/d_tau)
% beta=N*d_tau

d_tau=0.1;
N=round(max(beta)/d_tau)*ones(size(beta));
d_tau=beta./N

name="beta";
param = beta; %

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%stest_%s=%.15g.out',dossier,name,param(i)) ;
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('./%s%s config/transition.in %s=%.15g N_slices=%d output=%s', repertoire, code, name, param(i), N, output{i});
    disp(cmd)
%     system(cmd);

end

%% Analyse %%
%%%%%%%%%%%%%



H=cell(1,nsimul);
Hmoy=zeros(1,nsimul);

for i=1:nsimul
    data=load(output{i}+"_nrg.out");
    n_MCS=size(data,1);
    H{i}=data(:,1);
    Hmoy(i)=mean(H{i})/N(i);
end

%%
beta=beta(2:end);
Hmoy=Hmoy(2:end);


f=figure('Position',  [200, 200, 600, 350]);
plot(beta,Hmoy*1e-2,'+')
xlabel('$\beta_N$')
ylabel('$E \; {\rm [10^{-18} \; J]}$')
rm_space(gca);
print(f,'FIG/transition', '-depsc');

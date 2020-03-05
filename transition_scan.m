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

nsimul = 20; % number of simulations


beta=linspace(0.1, 30, nsimul);

d_tau=0.05;
N=round(beta/d_tau)
beta=N*d_tau


name="beta";
param = beta; %

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%s%s=%.15g.out',dossier,name,param(i)) ;
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('./%s%s config/transition.in %s=%.15g N_slices=%d output=%s', repertoire, code, name, param(i), N(i), output{i});
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
    Hmoy(i)=mean(H{i});
end

%%


figure
plot(beta,Hmoy,'+')
xlabel('$\beta$')
xlabel('$E$')

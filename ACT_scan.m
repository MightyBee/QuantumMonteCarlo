set(groot,'defaultAxesFontSize',16)
set(groot,'defaultAxesLabelFontSizeMultiplier',19/13)
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '';
code = 'main2';
dossier='simulations/ACT_scan/';

nsimul = 5; % number of simulations


N_slices=round(logspace(0.5, 2, nsimul));

name="N_slices";
param = N_slices; %


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%s%s=%.15g.out',dossier,name,param(i)) ;
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('./%s%s config/harmonic.in %s=%.15g output=%s', repertoire, code, name, param(i),output{i});
    disp(cmd)
    system(cmd);

end

%% Analyse %%
%%%%%%%%%%%%%



tMC_norm=cell(1,nsimul);
ACT_norm=cell(1,nsimul);
tau0=zeros(1,nsimul);
Neff=zeros(1,nsimul);

for i=1:nsimul

    data=load(output{i}+"_pos.out");
    n_part=1;
    n_MCS=size(data,1);
    n_slices=size(data,2)/n_part;
    A=zeros(n_MCS,n_slices,n_part);
    for k=1:n_part
        A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
    end
    Oi=mean(A(:,:,1),2);
    tMC=0:n_MCS-1;
    ACT=zeros(size(tMC));
    for k=1:n_MCS-1
        ACT(k)=1/(n_MCS-tMC(k)-1)*sum((Oi(1:n_MCS-tMC(k))-mean(Oi(1:n_MCS-tMC(k)))).*(Oi(tMC(k)+1:n_MCS)-mean(Oi(tMC(k)+1:n_MCS))));
        if ACT(k) <= 0
            iFinal=k-1;
            break;
        end

    end
    ACT0=ACT(1);
    ACT_norm{i}=ACT(2:iFinal)/ACT0;
    tMC_norm{i}=tMC(2:iFinal);

    tau0(i)=0.5+sum(ACT_norm{i})
    Neff(i)=n_MCS/(2*tau0(i))
end

%%
deltaT=1./N_slices

figure
loglog(1./N_slices,tau0,'+')

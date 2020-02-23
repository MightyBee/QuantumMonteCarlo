set(groot,'defaultAxesFontSize',16)
set(groot,'defaultAxesLabelFontSizeMultiplier',19/13)
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '../'; 
code = 'main2'; 
dossier='../simulations/ACT_scan/';

nsimul = 10; % number of simulations
 

N_slices=round(logspace(0.5, 2, nsimul));

name="N_slices";
param = N_slices; % 


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%s%s=%.15g.out',dossier,name,param(i)) ;
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('./%s%s ../config/harmonic.in %s=%.15g output=%s', repertoire, code, name, param(i),output{i});
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

% 
% t=cell(1,nsimul);
% h_moy=cell(1,nsimul);
% ht_moy=zeros(1,nsimul);
% sigma_h=zeros(1,nsimul);
% h_max=zeros(1,nsimul);
% P_h=cell(1,nsimul);
% s_ava=cell(1,nsimul);
% P_N=cell(1,nsimul);
% tc_moy=zeros(1,nsimul);
% sk=zeros(nsimul,4);
% 
% for i = 1:nsimul 
%     data = load(sprintf('%s_h.out', output{i}));
%     h=data;
%     t_fin=size(h,2);
%     t{i}=[1:t_fin]';
%     M=size(h,1);
%     h_moy{i}=zeros(1,t_fin);
%     for m=1:M
%         h_moy{i}=h_moy{i}+h(m,:);
%     end
%     h_moy{i}=h_moy{i}/M;
%     
%     data = load(sprintf('%s_t.out', output{i}));
%     t_c=data;
%     tc_moy(i)=mean(t_c);
%     
%     h_rec=h(1,t_c(1)+1:end);
%     ht_moy(i)=mean(h_rec);
%     sigma_h(i)=sqrt(mean(h_rec.^2)-ht_moy(i).^2);
%     
%     h_max(i)=max(h_rec);
%     P_h{i} = zeros(1,h_max(i));
%     for n = 1:h_max(i)
%         P_h{i}(n) = sum(h_rec==n);
%     end
%     P_h{i}=P_h{i}/length(h_rec);
%     tot=sum(P_h{i})
%     
%     data = load(sprintf('%s_bin.out', output{i}));
%     s_ava{i}=data(:,1);
%     P_N{i}=data(:,2);
%     
%     data = load(sprintf('%s_ava.out', output{i}));
%     ava=data(1,:);
%     for k = 1:4
%         sk(i,k)=mean(ava.^k);
%     end
% end
% 
% 
% 
% 
% %% Figures %%
% %%%%%%%%%%%%%
% k=1;
% fig(k)=figure('Position',[50,50,900,550]);
% clear p;
% for i=nsimul:-1:1
%     p(nsimul+1-i)=plot(t{i},h_moy{i},'DisplayName',sprintf('$L=%d$',L(i)));
%     hold on
% end
% hold off
% xlabel('$t$')
% ylabel('$\tilde{h}(t;L)$')
% ylim([0 1000])
% lgd=legend(p);
% set(lgd,'fontsize',16,'Location','northwest');
% print(fig(k),'figures/h0_moy_linlin', '-depsc');
% 
% k=k+1;
% fig(k)=figure('Position',[50,50,900,550]);
% clear p;
% for i=nsimul:-1:1
%     p(nsimul+1-i)=loglog(t{i}/L(i)^2,h_moy{i}/L(i),'DisplayName',sprintf('$L=%d$',L(i)));
%     hold on
% end
% hold off
% xlabel('$t/L^2$')
% ylabel('$\tilde{h}(t;L)/L$')
% xlim([2e-6 8e4])
% ylim([1e-3 3])
% lgd=legend(p);
% set(lgd,'fontsize',16,'Location','southeast');
% print(fig(k),'figures/h0_moy_collapse', '-depsc');
% 

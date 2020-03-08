set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultFigurePosition',  [200, 200, 450, 300]);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = '';
code = 'main3';
dossier='simulations/OH_scan/';

nsimul = 17; % number of simulations


xc2=(linspace(2.28, 3.03, nsimul));

name="xc2";
param = xc2; %


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%s%s=%.15g.out',dossier,name,param(i)) ;
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('./%s%s config/Hbonds.in %s=%.15g output=%s', repertoire, code, name, param(i),output{i});
    disp(cmd)
%     system(cmd);

end

%% Analyse %%
%%%%%%%%%%%%%



Heth=zeros(1,nsimul);
errHeth=zeros(1,nsimul);
Rmoy=zeros(1,nsimul);

for i=1:nsimul
    data=load(output{i}+"_e0.out");
    Heth(i)=data(1,1);
    errHeth(i)=data(1,2);
    data=load(output{i}+"_pos.out");
    n_part=3;
    n_MCS=size(data,1);
    n_slices=size(data,2)/n_part;
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
    sa=mean(B(:,2));
    B=B-sa;
    
    mean_x=zeros(1,3);
    std_x=zeros(1,3);
    for k=1:3
        mean_x(k)=mean(B(:,k));
        std_x(k)=std(B(:,k));
    end
    Rmoy(i)=mean_x(3)-mean_x(1)
    
%     figure('Position',[200,400,1520,680])
%     hold on;
%     for k=1:n_part
%         h(k)=histogram(B(:,k)-mean_x(1),-0.5:0.02:3,'Normalization','pdf');
% %         plot([mean_x(k) mean_x(k)]-mean_x(1),[0 25])
%     end
%     rm_space(gca);
end

%%

Rmoy=[2.0061 2.0556 2.1037 2.1535 2.1990 2.2495 2.2989 2.3477 2.3990 2.4538 2.5174 2.5867 2.6550 2.7381 2.8211 2.8986 2.9731]

errx=sqrt(std_x(1)^2+std_x(3)^2)*ones(1,nsimul);

figure
errorbar(Rmoy,Heth,errHeth,errHeth,errx,errx,'o')
xlabel('$R \; {\rm[\AA]}$');
ylabel('$E_0 \; {\rm[10^{-20} \, J]}$');
rm_space(gca);
%                             ????              ????  ????  ????                          ????
rmax_l(1,:)=     [1.00, 1.02, 1.04, 1.08, 1.10, 1.12, 1.14, 1.16, 1.04, 1.02, 1.00, 1.00, 0.98, 0.98, 0.98, 1.00, 0.98];
rmax_l(2,:)=Rmoy-[1.00, 1.02, 1.04, 1.08, 1.10, 1.12, 1.14, 1.16, 1.34, 1.42, 1.52, 1.60, 1.66, 1.74, 1.82, 1.90, 1.98];
rmax_r=rmax_l+0.02;
rmax=(rmax_l+rmax_r)/2;
% rmax=min(rmax,Rmoy-rmax)
rmax_lavie_l=max(rmax_l(1,:),rmax_l(2,:))-0.01;
rmax_lavie_r=min(rmax_r(1,:),rmax_r(2,:))+0.01;
err_lavie=(rmax_lavie_r-rmax_lavie_l)/2;
rmax_lavie=(rmax_lavie_r+rmax_lavie_l)/2;
erry=0.02*ones(1,nsimul);

figure
% errorbar(Rmoy,rmax(1,:),erry,erry,errx,errx,'o')
% hold on
% errorbar(Rmoy,rmax(2,:),erry,erry,errx,errx,'o')
% hold on
errorbar(Rmoy,rmax_lavie,err_lavie,err_lavie,errx,errx,'o')
hold on
xlabel('$R \; {\rm[\AA]}$');
ylabel('$r_m \; {\rm[\AA]}$');

Rstart = 2;
Rend = 3;

global R;
R = Rstart - (Rend - Rstart)/100;

x = linspace(Rstart, Rend);
y = x;

for i=1:100
    R = R + (Rend - Rstart)/100;
    y(i) = R/2 + fminbnd(@morse, -1, 0);
end;

plot(x, y, "Color", "b");
rm_space(gca);


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



data = load("output2_pos.out");

beta = 10;

hbar = 1.0545718;

R_reduced = [2.9 2.5 2.45 2.3];
R_list = [2.9 2.5 2.45 2.3 2.4 2.6 2.7 2.8 2.2 2.35 2.1];
E = 1.98630 * [-42 -47 -49 -57.27]; %e-20 J
% -83.4246  -93.3561  -97.3287 -113.7554

%n_thermal = 10 because we saw up to 6 in <x2>
% sometimes it was not symmetric
EH = [
-83.207152136852
-93.5683482384089
-97.0978470009506
-112.545840303505
-101.49572985321
-89.0594449137312
-86.0742398230457
-84.3727167751251
-125.760487408964
-106.716517821947
-140.478737208659
]; %e-20

EH_err = [
0.0703984067837869
0.0711027538630522
0.049224293870945
0.0431843331529173
0.0376806893087406
0.0692651293775684
0.07009193173918
0.0573431372840057
0.0432471126461167
0.0341888406718646
0.0526479453449164
]; %e-20

HH = [
-32.6771779116752
-41.5738060894126
-45.7008020177809
-59.7296699018952
-50.3073303690072
-37.0675626186829
-35.2785686503673
-32.3354113597087
-74.7406238432023
-54.0651065748491
-88.8698992298585
]; %e-20

HH_err = [
0.550409390908296
0.559774256573301
0.506940572585432
0.522706836480364
0.550059415265635
0.53302891595497
0.486308868651724
0.53768948377465
0.580446752846014
0.54481115983829
0.606178753479858
];


%-------------------------------------------


ED = [
-83.6931963902044
-94.1726641838847
-97.3551173194804
-112.661636616692
-101.553099624716
-89.4607823451549
-86.5428181503251
-84.7912319841338
-125.819801045429
-106.775596675297
-140.700154576665
]; %e-20

ED_err = [
0.0647844280833208
0.0627345853874779
0.0505084143983433
0.0375923031954617
0.0385336702907993
0.0644725141809373
0.0602306520383869
0.0639931847064842
0.0501195380143349
0.0402168345497726
0.0596400226703871
]; %e-20

HD = [
-31.6801007544756
-42.438295975053
-44.8731204757362
-61.5643728819758
-49.6712132805306
-37.5366404394425
-35.6678918159081
-33.0399398240751
-73.0447102593811
-54.4049936894798
-88.4331136684928
]; %e-20

HD_err = [
0.588330527444869
0.513685446695349
0.511976824487968
0.542924269109719
0.577234343138385
0.537944614304622
0.474967104925388
0.4982284640192
0.497077589398404
0.558348995127316
0.5408851880207
]; %e-20

%100 bins
rH_min_min = (-1/10000) * [4948 2611 2076 -12 1385 3197 3793 4364 0 918 0];
rH_min_max = (-1/10000) * [4746 2414 1878 127 1180 2984 3582 4136 0 679 0];

rD_min_min = (-1/10000) * [5000 2536 2297 0 1059 3296 3834 4591 0 920 0];
rD_min_max = (-1/10000) * [4800 2280 2076 0 856 3074 3635 4342 0 685 0];


rH_min = 0.5*(rH_min_max + rH_min_min);
rH_min_err  = 0.5*abs(rH_min_max - rH_min_min);

rD_min = 0.5*(rD_min_max + rD_min_min);
rD_min_err  = 0.5*abs(rD_min_max - rD_min_min);

figure('DefaultAxesFontSize', 14);


%{

plot(R_list, EH, '.', 'MarkerSize', 20, "Color", "r");
hold on;
plot(R_list, ED, '.', 'MarkerSize', 20, "Color", "g");
plot(R_reduced, E, '+', 'MarkerSize', 15, "Color", "k", 'LineWidth', 1);
%errorbar(R, EH, EH_err, 'LineStyle', 'none');

xlabel('$R$ [\AA]', 'Interpreter','latex', 'FontSize', 14);
ylabel('$E_0$ [$10^{-20}$ J]', 'Interpreter','latex','FontSize', 14);

%legend("Measured energies", "Reference energies", 'Interpreter','latex', 'FontSize', 14, 'Location', 'southeast');

ylim([-145 -80]);
yticks([-145:10:-80]);
xticks([2:0.1:3]);
xlim([2 3]);
%}





plot(R_list, R_list/2 + rH_min, '.', 'MarkerSize', 15, "Color", "r");
hold on;
errorbar(R_list, R_list/2 + rH_min, rH_min_err, 'LineStyle', 'none', "Color", "r");
plot(R_list, R_list/2 + rD_min, '.', 'MarkerSize', 15, "Color", "g");
errorbar(R_list, R_list/2 + rD_min, rD_min_err, 'LineStyle', 'none', "Color", "g");

xlabel('$R$ [\AA]', 'Interpreter','latex', 'FontSize', 14);
ylabel('$r_m$ [\AA]', 'Interpreter','latex','FontSize', 14);

%legend("Measured energies", "Reference energies", 'Interpreter','latex', 'FontSize', 14, 'Location', 'southeast');

ylim([0.9 1.25]);
yticks([0.9:0.05:1.25]);
xticks([2:0.1:3]);
xlim([2 3]);

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


%}






set(gcf, 'Position',  [200, 200, 450, 300]);

% magic code to remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

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

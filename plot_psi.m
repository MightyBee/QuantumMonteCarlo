data = load("output2_pos.out");

beta = 10;
hbar = 1.0545718;

xl = 1;
dx = 0.2;
yl = 2.5;
dy = 0.5;


%{
D = 83.402;
a = 2.2;
r0 = 0.96;
delta1 = 0.4*D;
b = 2.2;
R1 = 2*r0+1/a;
R = 2.9;
DELTA = delta1*exp(-b*(R-R1));

%----------------------------------
%Credits: https://www.chebfun.org/examples/ode-eig/DoubleWell.html

a = 1.0545718^2/(2*m)*100;

tic
x = chebfun('x');
L = chebop(-1,1);
%L.op = @(x,u) -a*diff(u,2) + (0.5*((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))+(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))) - sqrt(((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))-(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))))^2+4*DELTA^2))).*u;
L.bc = 0;
[EV,D] = eigs(L,1);
disp(diag(D)), toc
%------------------------------------
plot(EV(:,1)/sum(EV(:,1)), "b-", 'linewidth', 1.5);
%}


%0.5*((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))+(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))) - sqrt(((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))-(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))))^2+4*DELTA^2));

%(D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))

%(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0))))

%(D*(exp(-2*a*(x-r0))-2*exp(-a*(x-r0))));

figure('DefaultAxesFontSize', 14);
histogram(data(:),100,'FaceColor','#EDB120','Normalization','pdf');

%xlabel('$x \, [\AA]$', 'Interpreter','latex', 'FontSize', 16);
%ylabel('Probability density', 'Interpreter','latex','FontSize', 16);

%xticks([-xl:dx:xl]);
%xlim([-xl xl]);
%yticks([0:dy:yl]);
%ylim([0 yl]);


set(gcf, 'Position',  [200, 200, 600, 400]);

% magic code to remove white space
%{
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%}
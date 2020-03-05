data = load("output2_pot.out");
plot(2.2 * data(:,1), data(:,2)/83.402);
xlim([-2 2]);
%ylim([0 20]);

% 1.98630 e-23


%{
x = [-2:0.01:2];
D = 83.402;
a = 2.2;
r0 = 0.96;
delta1 = 0.4*D;
b = 2.2;
R1 = 2*r0+1/a;
R = 2.9;
DELTA = delta1*exp(-b*(R-R1));
y = (0.5*((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))+(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0)))) - sqrt(((D*(exp(-2*a*(R/2+x-r0))-2*exp(-a*(R/2+x-r0))))-(D*(exp(-2*a*(R/2-x-r0))-2*exp(-a*(R/2-x-r0))))).^2+4*DELTA^2)));
yy = 2*D*a*(exp(-a*(x-r0)) - exp(-2 *a*(x-r0)));

plot(x,y);
%hold on;
%plot(x,y);
%}
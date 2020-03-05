data = load("output2_pos.out");

figure('DefaultAxesFontSize', 12);

plot(mean(transpose(data.^2)));

m = 1;
w = 1;
beta = 0.5;
hbar = 1.0545718;

x = [-12:0.1:12];
%classic
y1 = sqrt(m*w^2*beta/(20*pi*hbar)) * exp(-m*w^2*beta/(20*hbar)* x.^2);
%quantum
y2 = sqrt(m*w/(10*pi*hbar)) * exp(-m*w*x.^2/(10*hbar));

plot(x, y1, "r-", "LineWidth", 1);
hold on;
plot(x, y2, "b-", "LineWidth", 1);
histogram(data(:),50,'FaceColor','#EDB120','Normalization','pdf');


xlabel('$x \, [\AA]$', 'Interpreter','latex', 'FontSize', 14);
ylabel('$|\psi(x)|^2$', 'Interpreter','latex','FontSize', 14);
legend("Thermal", "Quantum", "Simulation");

xticks([-12:2:12]);
xlim([-12 12]);
%yticks([0:0.02:0.2]);
%ylim([0 0.2]);


set(gcf, 'Position',  [200, 200, 500, 400]);


% magic code to remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];



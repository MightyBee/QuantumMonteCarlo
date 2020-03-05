n_masses = 10;
masses = logspace(-2, 2, n_masses);
times = [1 5 10];
%times = [1 5 10 15 20];
colors = ["r", "g", "b"];

%{
%results = zeros(n_masses,size(times,2));
masses = [1];
n_masses = 1;
results2 = zeros(n_masses,size(times,2));
for i=1:n_masses
    for j=1:size(times,2)
        %system("main2 configuration.in m1="+masses(i)+" beta="+times(j)+" N_slices="+(10*times(j)));
        system("main2 configuration.in m1=1 beta="+times(j)+" N_slices="+(10*times(j)));
        data = load("output2_pos.out");
        %results(i,j) = mean(std(transpose(data),1));
        results2(i,j) = mean(std(transpose(data),1));
    end;
end;
%}


figure('DefaultAxesFontSize', 14);

%loglog(times, results2, ".", 'MarkerSize', 20, "Color", "r");
%hold on;
%x = [1:0.1:20];
%loglog(x, 10^(-0.0626) * x.^0.5108, "--k", 'LineWidth', 1);
% 0.5108   -0.0626

x = [0.01:1:100];
for i=1:size(times,2)
    c = polyfit(log10(masses), log10(results(:,i).'), 1);
    disp(c);
    loglog(x, 10^c(2)*x.^c(1), "--k", "Color", colors(i), 'LineWidth', 1);
    hold on;
end;
for i=1:size(times,2)
    loglog(masses, results(:,i), ".", 'MarkerSize', 20, "Color", colors(i));
end;

legend("$\beta_N = 1, \,\,\,$ linear regression $y=-0.500x -0.060$", "$\beta_N = 5, \,\,\,$ linear regression $y=-0.498x + 0.303$", "$\beta_N = 10$, linear regression $y=-0.502x + 0.451$", 'Interpreter','latex', 'FontSize', 14);


%xlabel('$\beta$ [$10^{-15}$ s]', 'Interpreter','latex', 'FontSize', 14);
ylabel('$s$ [\AA]', 'Interpreter','latex','FontSize', 14);

%legend("Measured particle sizes", "Linear regression $y=0.511 x-0.063$", 'Interpreter','latex', 'FontSize', 14, 'Location', 'northwest');

xlabel('$m $ [$10^{-30}$ kg]', 'Interpreter','latex', 'FontSize', 14);

% for beta vs s
%xlim([0.95, 21]);

%xticks([0:1:10]);

set(gcf, 'Position',  [200, 200, 600, 400]);

% magic code to remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%data = load("output2_pot.out");
data = load("output2_pos.out");

line_width = 1;
colors = [
    [1 0 0];
    [0.9290, 0.6940, 0.1250];
    [0 1 0];
    [0 1 1];
    [0 0 1];
    [1 0 1];
    [0 0 0];
    ];

%show = mean(transpose(data)).^2 + std(transpose(data)).^2;
%plot(show);

figure('DefaultAxesFontSize', 12);

plot([0:1:9], data(1,:), '-', 'Color', colors(1,:), "LineWidth", line_width);
hold on;
grid on;
for i=2:size(data, 1)
    plot([0:1:9], data(i,:), '-', 'Color', colors(i,:), "LineWidth", line_width);
end;

xlabel('$k$', 'Interpreter','latex', 'FontSize', 14);
ylabel('$x_k \, [\AA]$', 'Interpreter','latex','FontSize', 14);

%xlim([0, 10]);
%xticks([0:1:10]);
%ylim([-1.5, 1.5]);

set(gcf, 'Position',  [200, 200, 600, 250]);

% magic code to remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
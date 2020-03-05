data = load("output2_pos.out");

speeds = diff(data(:,:).')/(0.1);

beta = 10;
hbar = 1.0545718;

figure('DefaultAxesFontSize', 14);
histogram(speeds,1000,'FaceColor','#EDB120','Normalization','pdf');

xlim([-1.2 1.2]);
%yticks([0:dy:yl]);
%ylim([0 yl]);

%xlabel('$x \, [\AA]$', 'Interpreter','latex', 'FontSize', 16);
%ylabel('Probability density', 'Interpreter','latex','FontSize', 16);

%xticks([-xl:dx:xl]);
%xlim([-xl xl]);
%yticks([0:dy:yl]);
%ylim([0 yl]);


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
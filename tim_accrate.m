set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/14)
set(0, 'DefaultLineLineWidth', 1.5);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% loading data %%


output='test_accrate';
accrate=load(output+"_rate.out");

%% figure %%
f=figure('Position',  [200, 200, 600, 350]);
plot((0:length(accrate)-1),accrate,'.')
xlabel('Metropolis Sweeps')
ylabel('Acceptance rate')
xlim([-2 100+2])
ylim([0 1])
% print(f,sprintf('FIG/x2_%s_%d',file,number), '-depsc');
% saveas(f,sprintf('FIG/x2_%s_%d.fig',file,number));


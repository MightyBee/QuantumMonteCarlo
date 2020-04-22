set(groot,'defaultAxesFontSize',13)
set(groot,'defaultAxesLabelFontSizeMultiplier',16/13)
set(0, 'DefaultFigurePosition',  [200, 200, 650, 400]);
set(0, 'DefaultLineMarkerSize', 11);
set(0, 'DefaultLineLineWidth', 1.2);
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%%
X=[-0.425992581468067;-0.712728711905450;-0.0602442913260656;-0.417749874673639;-0.567561250562530;-0.934579439252336;-0.920200449930611;-0.346467802463965;-0.496580803291655;0.0778114185860011];
X2=[0.183140859276149;0.741346260121131;0.807125917613432;0.278167771620860;-0.134579439252337;0.325693907296919;0.318767883956528;-0.0152999844592698;-0.154243234010262;-0.272076055876096];
n=length(X);
LM=X;
LM(5)=LM(5)-0.3;
GD=X+0.2;
BI=X;
BI(6:8)=BI(6:8)+0.4;
MI=-X;
CM=X-2*mean(X);
SW1=[X(2); X2(3:6); X(7)];
SW2=[X2(2); X(3:6); X2(7)];

yt=zeros(1,n+1);
ytl=cell(1,n+1);
for i=0:n
    yt(i+1)=i;
    ytl{i+1}=sprintf('%d',(mod(i,n)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=figure('Position',[900,500,950,800]);
% f1=figure('Position',[900,500,950,450]);

pos1 = [0.075 0.56 0.7/3 0.4];
subplot('Position',pos1)
plot([LM; LM(1)],[0:1:n],'--');
hold on
plot([X; X(1)],[0:1:n]);
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Local move')

pos2 = [0.175+0.7/3 0.56 0.7/3 0.4];
subplot('Position',pos2)
plot([GD; GD(1)],[0:1:n],'--');
hold on
plot([X; X(1)],[0:1:n]);
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Global displacement')

pos3 = [0.275+2*0.7/3 0.56 0.7/3 0.4];
subplot('Position',pos3)
plot([BI; BI(1)],[0:1:n],'--');
hold on
plot([X; X(1)],[0:1:n]);
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Bisection')
% print(f1,'FIG/moves1', '-depsc');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2=figure('Position',[900,500,950,450]);

pos1 = [0.075 0.06 0.7/3 0.4];
subplot('Position',pos1)
plot([MI; MI(1)],[0:1:n],'--');
hold on
plot([X; X(1)],[0:1:n]);
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Reflection')

pos2 = [0.175+0.7/3 0.06 0.7/3 0.4];
subplot('Position',pos2)
p(1)=plot([CM; CM(1)],[0:1:n],'--');
hold on
p(2)=plot([X; X(1)],[0:1:n]);
plot([mean(CM) mean(CM)],[0 n],':','Color',get(p(1),'Color'));
plot([mean(X) mean(X)],[0 n],':','Color',get(p(2),'Color'))
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Center of mass displacement')

pos3 = [0.275+2*0.7/3 0.06 0.7/3 0.4];
subplot('Position',pos3)
C=get(gca,'ColorOrder');
plot([X; X(1)],[0:1:n],'Color',C(2,:));
hold on
plot([X2; X2(1)],[0:1:n],'Color',C(3,:));
plot(SW1,[1:6],'--','Color',C(1,:));
plot(SW2,[1:6],'--','Color',C(5,:));
yticks(yt)
yticklabels(ytl)
xlim([-1 1])
xlabel('$x$')
ylabel('$\tau / \delta \tau$')
grid on
title('Swap')
% print(f2,'FIG/moves2', '-depsc');


print(f,'FIG/moves', '-depsc');

data = load("local_beta12_pos.out");

%{
norms = zeros(size(data,1),1);
%norms = data(:,1);

for i=1:size(data,1)
    %norms(i,1) = norm(data(i,:), 1);
    norms(i,1) = sqrt(data(i,1)^2 + data(i,2)^2);
end;

figure
histogram(norms,100,'Normalization','probability');
%}

%{
figure

t = plot([0:1:size(data,2)-1], data(1,:));
xlim([0 size(data,2)])
ylim([-10 10])

grid on

closw=waitforbuttonpress;

for i=2:size(data,1)
    pause(.0001)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'YData', data(i,:))
end

%}

%hbar	1.0545718 e-34
%kB		1.3806485 e-23

m = 1;
w = 1;
beta = 12;
hbar = 1.0545718;

x = [-10:0.1:10];
y1 = 10^(34/2-10) * sqrt(m*w^2*beta/(2*pi*hbar)) * exp(-m*w^2/(2*hbar)*10^14 * x.^2);
y2 = 0.5 * m * w * w * x.^2;
y3 = 10^(19/2-10) *sqrt(m*w/(pi*hbar)) * exp(-m*w^2*x.^2*0.1/hbar);
%%
figure
hold on;
% plot(x,y1, "r");
yyaxis left
plot(x,y3);
t=histogram(data(:),150,'FaceColor','#EDB120','Normalization','pdf');
xlim([-10 10])
xlabel('$x$')
ylabel('PDF')
yyaxis right
plot(x,y2);
ylabel('$V$')
%{
closw=waitforbuttonpress;
for i=2:1:size(data,1)
    pause(0.001)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'Data',data(i,:))
end
%}

%data = load("output2_pos.out");

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
ylim([-10 10])
xlim([0 size(data,2)])

ylabel('x')
xlabel('\tau / \delta \tau')

grid on

closw=waitforbuttonpress;

for i=2:size(data,1)
    pause(.1)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'YData', data(i,:))
end
%}

m = 10;
w = 1;
beta = 0.1;
T = 1/beta;

x = [-1:0.01:1];
y1 = sqrt(m*w*w/(2*pi*T)) * exp(-m*w*w/(2*T) * x.^2);
y2 = 0.5 * m * w * w * x.^2;

figure
hold on;
plot(x,y1, "r");
plot(x,y2, "b");
t=histogram(data(1,:),200,'Normalization','pdf');

closw=waitforbuttonpress;


for i=2:5:size(data,1)
    pause(0.01)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'Data',data(i,:))
end


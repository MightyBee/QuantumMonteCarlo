%data = load("output.out");

norms = zeros(size(data,1),1);
%norms = data(:,1);

for i=1:size(data,1)
    %norms(i,1) = norm(data(i,:), 1);
    norms(i,1) = sqrt(data(i,1)^2 + data(i,2)^2);
end;

figure
histogram(norms,100,'Normalization','probability');

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
    pause(.000001)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'YData', data(i,:))
end
%}



%{
figure
t=histogram(data(1,:),25,'Normalization','probability');
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:5:size(data,1)
    pause(.01)
    if ~ishandle(t)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(t,'Data',data(i,:))
end
%}

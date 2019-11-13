%{
data=load("output.out");
figure

t=plot([data(1,:) data(1,1)],[0:1:size(data,2)]);
xlim([-10 10])
xlabel('x')
ylabel('\tau / \delta \tau')
grid on
closw=waitforbuttonpress;


for i=2:size(data,1)
		pause(.01)
		if ~ishandle(t)
				break % Arrete l'animation si la fenetre est fermee
		end
		set(t,'XData',[data(i,:) data(i,1)])
end
%}

%%
n_part=2;
data=load("output.out");
%for k=1:n_part
%	A(:,:,k)=data(
%end
figure
t=histogram(data(1,:),91,'Normalization','probability');
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:size(data,1)
	pause(.01)
		if ~ishandle(t)
				break % Arrete l'animation si la fenetre est fermee
		end
		set(t,'Data',data(i,:))
end

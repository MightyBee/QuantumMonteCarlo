%% loading data %%
n_part=3;
data=load("output.out");
n_MCS=size(data,1);
n_slices=size(data,2)/n_part;
A=zeros(n_MCS,n_slices,n_part);
for k=1:n_part
	A(:,:,k)=data(:,(1+(k-1)*n_slices):k*n_slices);
end

%% imaginary path %%
%{
figure
k=3; % number of the particle that we want to analyse 
t=plot([A(1,:,k) A(1,1,k)],[0:1:size(A,2)]);
xlim([-10 10])
xlabel('x')
ylabel('\tau / \delta \tau')
grid on
closw=waitforbuttonpress;


for i=2:n_MCS
		pause(.01)
		if ~ishandle(t)
				break % Arrete l'animation si la fenetre est fermee
		end
		set(t,'XData',[A(i,:,k) A(i,1,k)])
end
%}

%% histogram %%
figure
hold on;
for k=1:n_part
    t(k)=histogram(A(1,:,k),91,'Normalization','probability');
end
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:n_MCS
	pause(.01)
		if ~ishandle(t(1))
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'Data',A(i,:,k))
        end
        
end



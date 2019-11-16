%% loading data %%
n_part=1;
data=load("output.out");
n_MCS=size(data,1);
n_slices=size(data,2)/(3*n_part);
A=zeros(n_MCS,n_slices,n_part,3);
for k=1:n_part
    for j=1:3
        A(:,:,k,j)=data(:,(j+3*(k-1)*n_slices):3:3*(k*n_slices-1)+j);
    end
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
    t(k)=histogram(sqrt(A(1,:,k,1).^2+A(1,:,k,2).^2+A(1,:,k,3).^2),91,'Normalization','probability');
    %t(k)=histogram(A(1,:,k,1),91,'Normalization','probability');
end
ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:n_MCS
	pause(.01)
		if ~ishandle(t(1))
				break % Arrete l'animation si la fenetre est fermee
        end
        for k=1:n_part
            set(t(k),'Data',sqrt(A(i,:,k,1).^2+A(i,:,k,2).^2+A(i,:,k,3).^2))
            %set(t(k),'Data',A(i,:,k,1))
        end
        
end

%% histogram %%
figure
hold on;
k=1;

for j=1:3
    t(j)=histogram(A(1,:,k,j),91,'Normalization','probability');
end

ylim([0 0.2])

closw=waitforbuttonpress;


for i=2:n_MCS
	pause(.01)
		if ~ishandle(t(1))
            break % Arrete l'animation si la fenetre est fermee
        end
        for j=1:3
            set(t(j),'Data',A(i,:,k,j))
        end
        
end



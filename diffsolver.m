%Credits: https://www.chebfun.org/examples/ode-eig/DoubleWell.html

figure;
hold on;

a = 1.0545718^2/2*100;

tic
x = chebfun('x');
L = chebop(-20,20);
L.op = @(x,u) -a*diff(u,2) + (6*(1-(x/6)^2)^2).*u;
L.bc = 0;
neigs = 1;
[EV,D] = eigs(L,neigs);
disp(diag(D)), toc


for j = 1:neigs
  d = D(j,j);
  plot(EV(:,j)/sum(EV(:,j)), 'LineWidth', 1.5, 'b-')
end
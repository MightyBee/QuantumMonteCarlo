data = load("output2_pos.out");

figure('DefaultAxesFontSize', 14);

%plot(mean(transpose(data.^2)));

m = 1;
w = 1;
L = 20;
V0 = 6;
x0 = 6;

beta = 10;
hbar = 1.0545718;

x = [-12:0.01:12];

xl = 12;
dx = 2;
yl = 0.22;
dy = 0.02;

%----------------------------------
%Credits: https://www.chebfun.org/examples/ode-eig/DoubleWell.html

a = 1.0545718^2/(2*m)*100;

tic
z = chebfun('z');
L = chebop(-25,25);
L.op = @(z,u) -a*diff(u,2) + (6*(1-(z/6)^2)^2).*u;
L.bc = 0;
[EV,D] = eigs(L,1);
disp(diag(D)), toc
%------------------------------------


%---classic
%double
fun = @(x) exp(-beta * V0 * (1 - (x/x0).^2).^2/(10 * hbar));
Z = integral(fun,-Inf,Inf);
y1 = exp(-beta * V0 * (1 - (x/x0).^2).^2/(10 * hbar))/Z;

%well
%y1 = probability_C(x, L);

%harmonic
%y1 = sqrt(m*w^2*beta/(20*pi*hbar)) * exp(-m*w^2*beta/(20*hbar)* x.^2);



%---quantum

%well
%y2 = probability_Q(x, L);

%harmonic
%y2 = sqrt(m*w/(10*pi*hbar)) * exp(-m*w*x.^2/(10*hbar));

plot(x, y1, "r-", "LineWidth", 1.5);
hold on;
%plot(x, y2, "b-", "LineWidth", 1.5);
plot(EV(:,1)/sum(EV(:,1)), "b-", 'linewidth', 1.5);
histogram(data(:),50,'FaceColor','#EDB120','Normalization','pdf');


xlabel('$x \, [\AA]$', 'Interpreter','latex', 'FontSize', 16);
ylabel('Probability density', 'Interpreter','latex','FontSize', 16);
%ylabel('$|\psi(x)|^2$', 'Interpreter','latex','FontSize', 14);
%legend("Classical mechanics", "Quantum", "Simulation");
%legend("Thermal", "Simulation");

dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
legend("$p_C(x)$", "$|\psi(x)|^2$", "PIMC", "$\beta_N = "+beta+"$", 'Interpreter','latex');

xticks([-xl:dx:xl]);
xlim([-xl xl]);
yticks([0:dy:yl]);
ylim([0 yl]);


set(gcf, 'Position',  [200, 200, 500, 400]);


%{
V0 = V0 * 1e-20;
x0 = x0 * 1e-10;

h = (8/x0^2)^(1/4);
c = sqrt(2 * V0/x0^4);

E = -h^8/(2^5 * c^2) + h^2/sqrt(2) - 2*c^2/h^4 - sqrt(2) * c^4 * 36/(8 * h^10);
E = E - 2^2 * h^2 * sqrt(h^6/(2*c^2))/(sqrt(pi) * 2^(1/4)) * exp(-h^6/(6*sqrt(2)*c^2));
%     ^^ +/-
%}

% magic code to remove white space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

function y = probability_Q(x, l)
    y = zeros(size(x));
    for i=1:size(y, 2)
        if(abs(x(i)) > l/2)
            y(i) = 0;
        else
            y(i) = 2/l * cos(pi * x(i)/l).^2;
        end
    end
end

function y = probability_C(x, l)
    y = zeros(size(x));
    for i=1:size(y, 2)
        if(abs(x(i)) > l/2)
            y(i) = 0;
        else
            y(i) = 1/l;
        end
    end
end



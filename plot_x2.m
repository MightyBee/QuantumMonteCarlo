data = load("output2_pos.out");

show = mean(transpose(data)).^2 + std(transpose(data)).^2;
plot(show);
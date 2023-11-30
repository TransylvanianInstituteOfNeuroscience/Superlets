slt20 = slt3(find(fois == 20), :);
slt40 = slt3(find(fois == 40), :);
slt60 = slt3(find(fois == 60), :);


figure;
hold on;
plot(time, slt20);
plot(time, slt40);
plot(time, slt60);
line([time(1), time(end)], [0.5, 0.5], 'color', 'black');


function f = plotcx(t, z)

f = figure;
hold on;
plot(t, real(z));
plot(t, imag(z));

return;
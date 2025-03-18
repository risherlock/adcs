clear all;


data = csvread('out_file.csv');
data = data(2:end, :);

time = data(:,1);
rx = data(:,2);
ry = data(:,3);
rz = data(:,4);

plot3(rx, ry, rz, 'r');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;


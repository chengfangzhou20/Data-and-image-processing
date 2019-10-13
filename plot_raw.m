%Plot data
clear all
clc
load Aged_8_amino.mat
for i = 3
figure(1)
hold on
plot(Ex{i},Sx{i},'-ob','LineWidth',2,'MarkerSize',10)
plot(Ey{i},Sy{i},'-or','LineWidth',2,'MarkerSize',10);
legend('Circumferential','Longitudianl')
end
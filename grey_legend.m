clc; close all; clear;

for idx = 1:7
greylevels = linspace(0,0.8,7); % for plotting

shade = [greylevels(idx) greylevels(idx) greylevels(idx)];

figure(1)
hold on
plot(1:7,idx*(1:7),'Color',shade,'LineStyle','-.','LineWidth',1.5)
end

figure(1)
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
legend('1','','','','','','7', 'FontName','Times New Roman', 'FontSize',18)
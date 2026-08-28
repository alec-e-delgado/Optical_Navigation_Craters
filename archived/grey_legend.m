clc; close all; clear;

n = 7;

% Create shades of green (dark -> light)
greens = [0 71 0; % dark
          0 117 0;
          0 163 0;
          0 209 0;
          0 255 0;
          138 255 138;
          184 255 184]/255; % light

% greylevels = linspace(0,0.8,7); % for grey shades

for idx = 1:7

shade = [greens(idx,1) greens(idx,2) greens(idx,3)];

figure(1)
hold on
plot(1:7,idx*(1:7),'Color',shade,'LineStyle','-.','LineWidth',1.5)
end

figure(1)
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
legend('1','','','','','','7', 'FontName','Times New Roman', 'FontSize',18)
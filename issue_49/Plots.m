close all;
HX_slices = 20;
WA1 = 5;
WA2 = 28;
WB1 = 46;
WBu = 46; 
t = 21601; 
tlim = 30;
font = 18; 
font2 = 16; 
Tmin = 0; 
Tmax = 320; 

cm=colormap(jet(t));
cm=flipud(cm);


time1 = (1 : 3036-2757+1) / 12; % There are 12, 5-min intervals in 1 hour 

% HELIUM IN
data1 = xlsread('XTData.xlsx', 1, 'C2757:C3036') + 273.15 * ones(3036-2757+1,1);

% He vs Time
figure 
hold on
plot((1 : t) * 1 / 3600, squeeze(T_data(end, 2, :, end)),...
        'k', 'LineWidth', 2)
plot(time1(1 : 12*6), data1(1 : 12*6), 'g', 'LineWidth', 2)    
set(gca,'fontsize',font2)
axis([0 6 0 305])
xlabel('Time (hr)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)
h_legend = legend('Model', 'Data', 'Location', 'Best');
set(h_legend,'FontSize',14);

% N2 vs Time
figure
plot((1 : t) * 1 / 3600, squeeze(q_data(end, :, end)),...
        'k','LineWidth',2)
set(gca,'fontsize',font2)
xlabel('Time (hr)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)

% Walls
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1680*(3/4), 350]);
    
subplot(1, 3, 1)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WA2 - WA1 - 1)*4.625*2.54*(WA2-WA1+1)/(WA2-WA1)^2, ...
        squeeze(T_wA_data(HX_slices/2, WA1 + 1 : WA2, :, 1)))
axis([0 4.625*2.54 Tmin Tmax])
set(gca,'fontsize', font2)
title('Cryopanel''s Shield, Left Side','FontSize', font)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)

subplot(1, 3, 2)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WBu - 1)*10*2.54*(WBu+1)/(WBu)^2,...
       squeeze(T_Bu_data(HX_slices/2, :, :)))
axis([0 10*2.54 Tmin Tmax])
set(gca,'fontsize',font2)
title('Cryopanel''s Bulk Head','FontSize', font)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)
    
subplot(1, 3, 3)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WA2 - WA1 - 1)*4.625*2.54*(WA2-WA1+1)/(WA2-WA1)^2, ...
        squeeze(T_wA_data(HX_slices/2, WA1 + 1 : WA2, :, 2)))
axis([0 4.625*2.54 Tmin Tmax])
set(gca,'fontsize',font2)
title('Cryopanel''s Shield, Right Side','FontSize', font)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)
    
set(gcf, 'Colormap', cm);
caxis([0 6])
hp4 = get(subplot(1,3,3),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  ...
    0.02  hp4(2)+hp4(3)*3.5], ...
    'Ticks',[0,1,2,3,4,5,6],...
    'TickLabels',{'0','1 hr','2 hr','3 hr','4 hr','5 hr','6 hr'}, 'FontSize', font2)

% Part 2
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1680*3/4, 350]);
subplot(1, 3, 1)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WA1 - 1)*2.54*(WA1 + 1)/WA1^2, ...
     squeeze(T_wA_data(HX_slices/2, 1 : WA1, :, 1)))
title('Cryopanel''s Wall, Left Side','FontSize', font)
axis([0 2.54 Tmin Tmax])
set(gca,'fontsize',font2)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)
    
subplot(1,3,2)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WB1 - 1)*10*2.54*(WB1+1)/WB1^2,...
        squeeze(T_wB_data(HX_slices/2, :, :)))
axis([0 10*2.54 Tmin Tmax])
set(gca,'fontsize',font2)
title('Cryopanel''s Wall, Bottom','FontSize', font)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font)
     
subplot(1,3,3)
set(gca, 'ColorOrder', cm);
hold all;
plot((0 : WA1 - 1)*2.54*(WA1 + 1)/WA1^2, ...
    squeeze(T_wA_data(HX_slices/2, 1 : WA1, :, 2)))
axis([0 2.54 Tmin Tmax])
set(gca,'fontsize',font2)
title('Cryopanel''s Wall, Right Side','FontSize', font)
xlabel('Length (cm)','FontSize', font)
ylabel('Temperature (K)','FontSize', font) 
    
set(gcf, 'Colormap', cm);
caxis([0 6])
hp4 = get(subplot(1,3,3),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  ...
    0.02  hp4(2)+hp4(3)*3.5], ...
    'Ticks',[0,1,2,3,4,5,6],...
    'TickLabels',{'0','1 hr','2 hr','3 hr','4 hr','5 hr','6 hr'}, 'FontSize', font2)

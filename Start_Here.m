

clc;
clear;
Text = sprintf('THS system initiated');
disp(Text);
pause(1)
Text = sprintf('Waiting for user...');
disp(Text);

Font= 'calibri';
FontSize = 8;
% set(0, 'DefaultUicontrolFontName', Font);
% set(0, 'DefaultUicontrolFontSize', FontSize + 4);
% set(0, 'DefaultUicontrolUnits','normalized');
% set(0, 'DefaultUicontrolHandleVisibility','off');
% set(0, 'DefaultUicontrolBackgroundColor','white');
% set(0, 'DefaultUipanelFontName', Font);
% set(0, 'DefaultUipanelHandleVisibility', 'off');
% set(0, 'DefaultUipanelFontSize', FontSize);
set(0, 'DefaultAxesFontName', Font);
set(0, 'DefaultAxesFontSize', FontSize);
% set(0, 'DefaultUitableUnits','normalized');
% set(0, 'DefaultUitableFontName', Font);
% set(0, 'DefaultUitableFontSize', FontSize-2);
set(0, 'DefaulttextFontName',Font)
set(0, 'DefaulttextFontSize',FontSize)

addpath(genpath([cd '\Model']));
addpath(genpath([cd '\DrivingCycle']));

% Driving Cycle Interface %
figure('Visible','on',...
    'Position',[450,150,450,300],...
    'Name','THS Optimal Control',...
    'NumberTitle','off',...
    'Color', 'white' ,...
    'MenuBar','none');
S = uicontrol('parent', gcf,'Style', 'listbox','units','normal',...
    'String', 'US06-HWY|EUDC|NEDC|FTP-75|HWFET|ECE-15|FTP-HWY|UDDS|10-Mode|15-Mode|HWFET-MTN|WLTC1|WLTC2|WLTC3',...
    'Position', [0.05 0.1 0.4 0.8]);
set(S,'Value',12);

D = uipanel('parent', gcf,...
    'Title','Tesla',...
    'Position',[0.5,0.32,0.45,0.6],...
    'BackgroundColor','white',...
    'Borderwidth',1);
axes('parent',D,...
    'Position',[0 0.05 0.98 0.9]);
[X,map] = imread('Model\Tesla.png','png');
imshow(X,map);

E = uipanel('parent', gcf,...
    'Title','Actions',...
    'Position',[0.5,0.1,0.45,0.2],...
    'BackgroundColor','white',...
    'Borderwidth',1);

F = uicontrol('parent',E,...
    'Style', 'pushbutton',...
    'String', '',...
    'Units','Normal',...
    'Position', [0.05 0.05 0.9 0.92],...
    'String','Run Simulation ',...
    'callback',@(src,event)runner(S,src,event));
%     'CData',imread('Model\play.png','png'),...


clear;





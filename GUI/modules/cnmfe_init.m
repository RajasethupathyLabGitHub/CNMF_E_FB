
panel_init = uipanel(...
'Parent',main_fig,...
'FontUnits',get(0,'defaultuipanelFontUnits'),...
'Units','characters',...
'Title','Initialization',...
'Position',[5 1 67 11],...
'BackgroundColor',get(0,'defaultuipanelBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuipanelBackgroundColorMode'),...
'ParentMode','manual',...
'Tag','data_panel',...
'FontSize',18,...
'FontWeight','bold');

push_init = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String','Initialize A & C',...
'Style','pushbutton',...
'Position',[18 7 30 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',16, 'fontweight', 'bold', ...
'callback', 'push_init_callback;');

% debug mode 
debug_on = true; 
save_avi = false; 

check_debug = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String','debug mode',...
'Style','checkbox',...
'value', 1, ...
'Position',[2 5 20 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14, 'fontweight', 'bold', ...
'callback', 'debug_on = get(check_debug, ''value'');');

check_saveavi = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String','save AVI',...
'Style','checkbox',...
'value', 0, ...
'Position',[23, 5, 20, 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14, 'fontweight', 'bold',...
'callback', ['save_avi=get(check_saveavi, ''value''); ',...
'debug_on=or(get(check_saveavi, ''value''), get(check_debug, ''value''));', ...
'set(check_debug, ''value'', debug_on);']);

%thresholds for detecting seed pixels 
text_mincorr = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 'min corr.',...
'Style','pushbutton',...
'Position',[2 3 15 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

edit_mincorr = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', num2str(default_options.min_corr'),...
'Style','edit',...
'Position',[17.5 3 10 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14,...
'callback', 'neuron.options.min_corr= str2double(get(edit_mincorr, ''string'')));');

text_minpnr = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 'min PNR.',...
'Style','pushbutton',...
'Position',[30 3 15 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

edit_minpnr = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', num2str(default_options.min_pnr'),...
'Style','edit',...
'Position',[45.5 3 10 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14,...
'callback', 'neuron.options.min_corr= str2double(get(edit_mincorr, ''string'')));');

%% 
text_patch = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', '# patch',...
'Style','pushbutton',...
'Position',[2 1 15 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

edit_nr = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 2,...
'Style','edit',...
'Position',[17.5 1 3.5 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

text_times = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 'X',...
'Style','text',...
'Position',[21 1 3 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

edit_nc = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 2,...
'Style','edit',...
'Position',[24 1 3.5 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);


text_patch = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 'max # neuron',...
'Style','pushbutton',...
'Position',[30 1 20 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);

edit_K = uicontrol(...
'Parent',panel_init,...
'FontUnits',get(0,'defaultuicontrolFontUnits'),...
'Units','characters',...
'String', 500,...
'Style','edit',...
'Position',[50 1 5 1.75],...
'BackgroundColor',get(0,'defaultuicontrolBackgroundColor'),...
'BackgroundColorMode',get(0,'defaultuicontrolBackgroundColorMode'),...
'Callback','%automatic',...
'Children',[],...
'ParentMode','manual',...
'Tag','data_dir',...
'FontSize',14);
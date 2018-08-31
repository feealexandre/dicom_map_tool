function [selection,value] = d_filedlg(varargin)
%DCE_FILEDLG IS A SHAMELESS RIP-OFF OF MATLAB's LISTDLG
% -- FONT IS DIFFERENT
% -- SOME FUNCTIONALITY IS LOST
% -- SOME FINCTIONALITY ADDED!!
%
%  See also LISTDLG

% Additions:
% 7/2013, Mikko Nissi
% - Changed promptstring fontsize to 14 BOLD to pull more attention
% - RESIZING!!!
%  Mikko Nissi, 2013

% RIP-OFF by Mikko Nissi, mikko.nissi@iki.fi

error(nargchk(1,inf,nargin))

figname = '';
smode = 2;   % (multiple)
promptstring = {};
liststring = [];
listsize = [240 550];
initialvalue = [];
okstring = 'OK';
cancelstring = 'Cancel';
fus = 8;
ffs = 8;
uh = 22;

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error('MATLAB:listdlg:InvalidArgument', 'Arguments to LISTDLG must come param/value in pairs.')
end
for i=1:2:length(varargin)
    switch lower(varargin{i})
     case 'name'
      figname = varargin{i+1};
     case 'promptstring'
      promptstring = varargin{i+1};
     case 'selectionmode'
      switch lower(varargin{i+1})
       case 'single'
        smode = 1;
       case 'multiple'
        smode = 2;
      end
     case 'listsize'
      listsize = varargin{i+1};
     case 'liststring'
      liststring = varargin{i+1};
     case 'initialvalue'
      initialvalue = varargin{i+1};
     case 'uh'
      uh = varargin{i+1};
     case 'fus'
      fus = varargin{i+1};
     case 'ffs'
      ffs = varargin{i+1};
     case 'okstring'
      okstring = varargin{i+1};
     case 'cancelstring'
      cancelstring = varargin{i+1};
     otherwise
      error('MATLAB:listdlg:UnknownParameter', ['Unknown parameter name passed to LISTDLG.  Name was ' varargin{i}])
    end
end

if ischar(promptstring)
    promptstring = cellstr(promptstring); 
end

if isempty(initialvalue)
    initialvalue = 1;
end

if isempty(liststring)
    error('MATLAB:listdlg:NeedParameter', 'ListString parameter is required.')
end

ex = get(0,'defaultuicontrolfontsize')*1.7;  % height extent per line of uicontrol text (approx)

fp = get(0,'defaultfigureposition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'defaultUicontrolBackgroundColor') ...
    'resize'                 'on' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'visible'                'off' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'ResizeFcn',             @l_ResizeFcn,...
    'closerequestfcn'        'delete(gcbf)' ...
            };
   % 'windowstyle'            'modal' ...

liststring=cellstr(liststring);

fig = figure(fig_props{:});

if length(promptstring)>0
    prompt_text = uicontrol('style','text','string',promptstring,...
        'horizontalalignment','left',...
        'fontsize',16,...
        'fontweight','bold',...
        'tag','promptstring',...
        'position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring)) ...
        listsize(1) ex*length(promptstring)]); %#ok
end

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;

listbox = uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
                    'string',liststring,...
                    'backgroundcolor','w',...
                    'max',smode,...
                    'fontname','courier','fontsize',10,...
                    'tag','listbox',...
                    'value',initialvalue, ...
                    'callback', {@doListboxClick});

ok_btn = uicontrol('style','pushbutton',...
                   'string',okstring,...
                   'position',[ffs+fus ffs+fus btn_wid uh],...
                   'callback',{@doOK,listbox});

cancel_btn = uicontrol('style','pushbutton',...
                       'string',cancelstring,...
                       'position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
                       'callback',{@doCancel,listbox});

if smode == 2
    selectall_btn = uicontrol('style','pushbutton',...
                              'string','Select all',...
                              'position',[ffs+fus 4*fus+ffs+uh listsize(1) uh],...
                              'tag','selectall_btn',...
                              'callback',{@doSelectAll, listbox});

    if length(initialvalue) == length(liststring)
        set(selectall_btn,'enable','off')
    end
    set(listbox,'callback',{@doListboxClick, selectall_btn})
end

set([fig, ok_btn, cancel_btn, listbox], 'keypressfcn', {@doKeypress, listbox});

% should be original version for this to work....:
%set(fig,'position',getnicedialoglocation(fp, get(fig,'Units')));
% Make ok_btn the default button.
%setdefaultbutton(fig, ok_btn);

% make sure we are on screen
movegui(fig)
set(fig, 'visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(listbox);
    uiwait(fig);
catch
    if ishandle(fig)
        delete(fig)
    end
end

if isappdata(0,'ListDialogAppData__')
    ad = getappdata(0,'ListDialogAppData__');
    selection = ad.selection;
    value = ad.value;
    rmappdata(0,'ListDialogAppData__')
else
    % figure was deleted
    selection = [];
    value = 0;
end

%% Resize function
function l_ResizeFcn(h,evd)
try
    screen=get(0,'screensize');
    H.FIG=gcbf;
    H.ListBox=findobj(get(H.FIG,'child'),'tag','listbox');
    H.PromptString=findobj(get(H.FIG,'child'),'tag','promptstring');
    
    fig_pos=get(H.FIG,'position');
            
    % check if minimum position is requirement is met
    if fig_pos(3)<240 fig_pos(3)=240; end
    if fig_pos(4)<300 fig_pos(4)=300; end
            
    if fig_pos(2)+fig_pos(4)>screen(4)-43
        % top of the figure is outside screen.. bring back.
        fig_pos(2)=screen(4)-fig_pos(4)-43;
    end
    set(H.FIG,'position',fig_pos);
    
    % save new position to preferences..
    Dat.HSiz=fig_pos(3)/screen(3);
    Dat.VSiz=fig_pos(4)/screen(4);
    Dat.HPos=fig_pos(1)/screen(3);
    Dat.VPos=fig_pos(2)/screen(4);
    
    setpref('d_FILEDLG','HSiz',Dat.HSiz);
    setpref('d_FILEDLG','VSiz',Dat.VSiz);
    setpref('d_FILEDLG','HPos',Dat.HPos);
    setpref('d_FILEDLG','VPos',Dat.VPos);
            
    figw=round(Dat.HSiz*screen(3)); % set screen to whatever it was
    figh=round(Dat.VSiz*screen(4)); % set screen to whatever it was
    dx=10;
    dy=10;
    
    % Panels & boxes & others..
    Dat.Pos.ListBox=[dx 5*dy fig_pos(3)-2*dx fig_pos(4)-8*dy];
    Dat.Pos.PromptString=[dx fig_pos(4)-3*dy fig_pos(3)-2*dx 3*dy];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resizes
    % Background panel for information and adjusts
    set(H.ListBox,'position',Dat.Pos.ListBox);
    set(H.PromptString,'position',Dat.Pos.PromptString);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
catch
    aedes_errordump(lasterror);
end
 % l_ResizeFcn




%% figure, OK and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
 case 'escape'
  doCancel([],[],listbox);
end

%% OK callback
function doOK(ok_btn, evd, listbox) %#ok
if (~isappdata(0, 'ListDialogAppData__'))
    ad.value = 1;
    ad.selection = get(listbox,'value');
    setappdata(0,'ListDialogAppData__',ad);
    delete(gcbf);
end

%% Cancel callback
function doCancel(cancel_btn, evd, listbox) %#ok
ad.value = 0;
ad.selection = [];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox) %#ok
set(selectall_btn,'enable','off')
set(listbox,'value',1:length(get(listbox,'string')));

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);
elseif nargin == 3
    if length(get(listbox,'string'))==length(get(listbox,'value'))
        set(selectall_btn,'enable','off')
    else
        set(selectall_btn,'enable','on')
    end
end

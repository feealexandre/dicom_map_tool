function dce_tablet(data,ttext)
%DCE_TABLET shows input matrix in a table
%    DCE_TABLET(DATA) shows 2D data in DATA. 
%    DATA can be a regular matrix or cell matrix
%    Displaying of numeric data can be adjusted from 
%    the GUI, others remain untreated.

% Changes
% 11/2013, Mikko Nissi
%   - added title to the inputs
%
% Mikko Nissi 2011, <mikko.nissi@iki.fi>

% Public variables
H=[];   
Dat=[];



% check input
if nargin<1
    error('No input!');
elseif nargin<2
    if iscell(data)
        Dat.Input=data;
    elseif ndims(data)<3
        % Make it to be a cell!
        Dat.Input=l_Cellize(data);
    else
        error('Error or cannot display more than 2D of data');
    end
    % Draw & update GUI
    Dat.ttext='Data tablet';
    l_DrawGUI;
    l_Update;
    
elseif nargin<3
    if isempty(ttext)
        Dat.ttext='Data tablet';
    else
        Dat.ttext=ttext;
    end
    if iscell(data)
        Dat.Input=data;
    elseif ndims(data)<3
        % Make it to be a cell!
        Dat.Input=l_Cellize(data);
    else
        error('Error or cannot display more than 2D of data');
    end
    
    % Draw GUI
    l_DrawGUI;
    l_Update;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Draw GUI objects
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_DrawGUI
  % Position figure slightly right of the center of the screen
  screen=get(0,'screensize');
  DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
  
  fig_w = 600;
  fig_h = 400;
  fig_pos = [min(screen(3)/2,screen(3)-fig_w) screen(4)/2-fig_h/2 fig_w fig_h];
  Dat.Position.Old=fig_pos;
  
  %% Main Figure ----------------------------
  H.MAINFIG = figure('Units','Pixel', ...
    'position',fig_pos,...
    'Name',Dat.ttext, ...
    'Numbertitle','off', ...
    'Tag','dce_tablet', ...
    'Color',DefaultColor, ...
    'Toolbar','none', ...
    'Menubar','none', ...
    'DoubleBuffer','on', ...
    'DockControls','off',...
    'renderer','painters',...
    'KeyPressFcn','',...
    'CloseRequestFcn',@l_Quit,...
    'Handlevisibility','off');
  
  H.CONTROLBOX=uipanel('parent',H.MAINFIG,...
    'units','pixel',...
    'position',[0 0 fig_w 50]);
  
  
  H.DECSEPARATOR_TX = uicontrol('parent',H.CONTROLBOX,...
    'units','pixel',...
    'position',[10 30 130 15],...
    'style','text',...
    'string','Decimal separator:',...
    'horizontalalign','left',...
    'backgroundcolor',DefaultColor);

  tmp = get(H.DECSEPARATOR_TX,'position');
	val=1;

  H.DECSEPARATOR = uicontrol('parent',H.CONTROLBOX,...
    'units','pixel',...
    'position',[110 33 100 15],...
    'style','popup',...
    'string',{'. (point)',', (comma)'},...
    'value',val,...
    'backgroundcolor','w',...
    'callback',@l_Update);
  
  H.NUMDEC_TX =  uicontrol('parent',H.CONTROLBOX,...
    'units','pixel',...
    'position',[10 10 100 15],...
    'style','text',...
    'string','Number of decimals:',...
    'horizontalalign','left',...
    'backgroundcolor',DefaultColor);

  val = '4';
  H.NUMDEC = uicontrol('parent',H.CONTROLBOX,...
    'units','pixel',...
    'position',[110 10 100 19],...
    'style','edit',...
    'string',val,...
    'backgroundcolor','w',...
    'callback',@l_Update);
  
  H.RESTABLE = uitable('Parent',H.MAINFIG,'position',[0 50 fig_w fig_h-51]);
  
  % Check Matlab version since uitable properties have changed in R2008a
  Dat.MatlabVersion = l_getmatlabversion;
  if Dat.MatlabVersion>=7.06
    set(H.RESTABLE,'Enable','on','visible','off');
  else
    set(H.RESTABLE,'Editable',false,'visible',false)
  end
  
  H.CLOSE=uicontrol('parent',H.CONTROLBOX,...
    'units','pixel',...
    'position',[fig_w-100 12 80 25],...
    'string','Close',...
    'callback',@l_Quit);
    
  % Set resize function
  set(H.MAINFIG,'resizefcn',@l_Resize)
  
end % function l_DrawGUI

%%
% Get matlabversion:
function version_number=l_getmatlabversion
  version_str = version;
  isImProcToolbox = false;
  ind = find(version_str=='.');
  major_ver = str2double(version_str(1:ind(1)-1));
  minor_ver = str2double(version_str(ind(1)+1:ind(2)-1));
  version_number = str2double(sprintf('%d.%02d',major_ver,minor_ver));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Resize window
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Resize(h,evd)
  fig_pos = get(H.MAINFIG,'position');
  if fig_pos(4)<250 || fig_pos(3)<320
    set(H.MAINFIG,'position',Dat.Position.Old)
    return
  end
  tmp=get(H.CONTROLBOX,'position');
  set(H.CONTROLBOX,'position',[0 0 fig_pos(3) 50]);
  controlpos=get(H.CONTROLBOX,'position');
  set(H.RESTABLE,'position',[0 50 fig_pos(3) fig_pos(4)-51]);
  set(H.CLOSE,'position',[fig_pos(3)-100 12 80 25]);  
  Dat.Position.Old=fig_pos;    
  l_Update([],[]);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Update table
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Update(h,evd)
  
  if isempty(Dat.Input)
    return
  end
  
  decsepval = get(H.DECSEPARATOR,'value');
  if decsepval==1
    decsep = '.';
  else
    decsep = ',';
  end
  numdec = get(H.NUMDEC,'string');
  numdec = str2num(numdec);
  numdec = floor(numdec);
  if isempty(numdec) || ~isreal(numdec) || numdec<0
    errordlg('The "Number of decimals" value must be a positive integer!',...
             'Error in "Number of decimals" field','modal')
    return
  end
  
  % Modify / check input cell
  grr=l_WorkInput(Dat.Input,decsep,numdec);
    
  
  Tab_h=size(grr,2);
  Tab_v=size(grr,1);
  fig_pos = get(H.MAINFIG,'position');
  
  % Adjust column width to something possibly reasonable..
  ColW = max(90,floor((fig_pos(3)-55) / Tab_h));
  
  
  % Clear table
  if Dat.MatlabVersion>=7.06
    set(H.RESTABLE,'data',[])
    set(H.RESTABLE,'Data',grr,...
      'visible','on',...
      'ColumnWidth',{ColW},... % repmat({90},1,Tab_h),...
      'ColumnFormat',repmat({'char'},1,Tab_h),...
      'ColumnName',[])
  else
    set(H.RESTABLE,'NumRows',2,...
      'NumColumns',2);
    drawnow
    drawnow
    tmp={'','';'',''};
    set(H.RESTABLE,'Data',tmp)
    drawnow
    drawnow
	
    % Update table
    set(H.RESTABLE,'visible',true)
    drawnow;
    drawnow;
    set(H.RESTABLE,'NumRows',Tab_v,'NumColumns',Tab_h);
    drawnow;
    drawnow;
    set(H.RESTABLE,'Data',grr);
    drawnow
    drawnow
  end
  
end % function l_Update


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adjust input data as requested!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=l_WorkInput(dat,decsep,numdec)
  f={''};
  if ~isempty(dat)
    for mm=1:size(dat,1)
      for nn=1:size(dat,2)
        tmp=dat{mm,nn};       
        if isnumeric(dat{mm,nn})
          if numdec>0       
            signd=sign(dat{mm,nn});
            floord=floor(abs(dat{mm,nn}));
            decs=(abs(dat{mm,nn})-floor(abs(dat{mm,nn})))*10^(numdec);
            f(mm+1,nn)={sprintf(['%1.0f' decsep '%0' num2str(numdec) '.0f'],signd*floord,decs)};
          else
            f(mm+1,nn)={sprintf('%1.0f',dat{mm,nn})};
          end
        else
          f(mm+1,nn)=dat(mm,nn);
        end        
      end
    end
  end
  
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Turn into cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=l_Cellize(data)
  f={''};
  if ~isempty(data)
    for mm=1:size(data,1)
      for nn=1:size(data,2)
        f(mm,nn)={data(mm,nn)};
      end
    end
  end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Quit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_Quit(h,evd)  
  % Clear variables
  h=H.MAINFIG;
  clear H Dat
  
  % Close window
  delete(h)
end

end

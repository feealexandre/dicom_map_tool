function f=d_loader(varargin)
%D_LOADER is an standalone / AEDES loader front-end for dicom files

% 11/18/2013, Mikko Nissi
%    - Fixed a bug related to opening single 4-D dataset with load_mMRI: if
%      only one file was selected, sorting by fip angle would cause most of
%      the data to be lost..
%      
% 7/30/2013, Mikko Nissi
%    - Fixed another bug from file loading and caching.. now checks whether
%      the files are dicom in addition to other checks. This should signifi
%      cantly aid keeping the cache current!
%    - Fixed a bug related to d_loader output, whether run with/without
%      output assigned
%
% 7/30/2013, Mikko Nissi
%    - Changed input to varargin to allow controlling running of the
%      function. The following modes implemented currently:
%          - 'openmode',[1,2,3] --> go directly to one of the existing
%          loaders and do not display the main figure
%          - 'title', STRING -> use STRING as the heading for file list!
%
% 7/25-27/2013, Mikko Nissi
%    - Modified the loading portion to allow loading of
%          - multiple 4-D datasets, as long as they have the same dimensions
%          - "dissimilar structures", i.e. multiple datasets with different
%            properties, again as long as the dimensions are the same
%          --> These modifications should allow loading data such as
%             1) MT/T1sat measurements with +Z and -Z acquired separetely
%             2) RAFF measurements similarly +Z/-Z separated
%             3) phase/magnitude data to be loaded
%
%    - Added possibility for function-use and returning data - in that
%      case, no Aedes is started!
%
% 7/2013, Mikko Nissi
%    - corrected a bug realated to loading an uncached directory
%
% 5/2013, Mikko Nissi
%    - Clean up of code.. Reduced number of loaders as they were redundant
%
% 3/2012, Mikko Nissi
%    - Due to gerenal usefulness, changed name to D_LOADER
%
% 3/2012, Mikko Nissi
%    - Several bug fixes and updates
%    - Most notably cache for directories added (significantly faster loading times after first open)
%      and reduced the amount of data saved in the index file (much faster
%      loading times)
%
% Mikko Nissi, 2011
%    - Updates, most of the program anyways..
%
% (c) Mikko Nissi, 2011 <mikko.nissi@iki.fi>


%% Public variables
H=[];   % handles to everything
Dat=[]; % data
Dat.d_indexfile='d_indexfile.mat'; % filename containing dicomindex
Dat.openmode=0; % default to regular run
Dat.title='Select file(s)';  % default promptstring for file choice list
nowait=0;   % dummy non-clearing variable to control waitfor if/not in openmode

% parse inputs
for tmpii=1:2:length(varargin) % scan every other input for tags
    switch varargin{tmpii}
        case 'openmode'
            Dat.openmode=varargin{tmpii+1};
            nowait=1;
        case 'title'
            Dat.title=varargin{tmpii+1};
    end
    % all should be set..
end

if nargout>0
    f=[];  % default to nothing
    Dat.To_Aedes=0;
else
    Dat.To_Aedes=1;
end


% Detect Matlab version
[Dat.MatlabVersion,Dat.isImageProc] = aedes_getmatlabversion;
% Check some preferences
l_Check_prefs;
% Draw front end GUI and clear all
H=l_draw_gui;
set(H.FIG,'ResizeFcn',@l_ResizeFcn);

l_ClearData;

if Dat.openmode>0
    set(H.FIG,'visible','off');
    switch Dat.openmode
        case 1
            l_Open_DATA(H.BTN_MRI,[],[]);
        case 2
            l_Open_DATA(H.BTN_mMRI,[],[]);
        case 3
            l_Open_DATA(H.BTN_PET,[],[]);            
    end
end

% Wait for quit if nargout
if nargout>0 && ~nowait
    waitfor(H.FIG);
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Debug function - paste debug information into workspace...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l_debug(h,evd)
        assignin('base','param',varargin);
        assignin('base','H',H);
        assignin('base','Dat',Dat);
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DRAW GUI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function H=l_draw_gui
        
        debug=true;
        
        % Determine paths
        tmp=which('d_loader.m');
        [fp,fn,fe]=fileparts(tmp);
        Dat.dFolder=[fp,filesep];
        
        % Load default font and colors
        if isunix
            DefaultColor = [239 235 222]/255;
        else
            DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
        end
        GD=aedes_gui_defaults;
        GD.col.frame = DefaultColor;
        
        % Figure Size Definitions ----------------
        screen=get(0,'ScreenSize'); % monitor resolution
        
        figw=round(Dat.HSiz*screen(3)); % set screen to whatever it was
        figh=round(Dat.VSiz*screen(4)); % set screen to whatever it was
        dx=10;
        dy=10;
        editw=45;
        edith=20;
        textw=120;
        axh=figh/2-2*dy; % axis height;
        axw=figw/2-2*dx; % axis width;
        Dat.Eps=1/axw; % "small nonzero unit"
        loadh=figh-8*dy;
        loadw=figw-2*dx;
        btnh=loadh/3-dy;  % 3 buttons..
        btnw=loadw-1.5*dx;
        
        
        
        % Panels & boxes & others..
        Dat.Pos.LoadPanel=[dx 7*dy loadw loadh];
        Dat.Pos.LoadMRI  =[dx/2 loadh-1*btnh-2*dy btnw btnh];
        Dat.Pos.LoadmMRI  =[dx/2 loadh-2*btnh-2*dy btnw btnh];
        Dat.Pos.LoadPET  =[dx/2 loadh-3*btnh-2*dy btnw btnh];
        Dat.Pos.Quit     =[dx*8 dy btnw-10*dx 5*dy];
        
        % Calculate gui position on the screen
        fig_pos = [round(Dat.HPos*screen(3)) round(Dat.VPos*screen(4)) figw figh];
        H.ORIG_FIG_POS = fig_pos;
        
        % Check if other windows exist
        figs=findall(0,'tag','d_loader_fig');
        if ~isempty(figs)
            tmp_pos=get(figs(1),'position');
            fig_pos(1)=tmp_pos(1)+15;
            fig_pos(2)=tmp_pos(2)-15;
        end
        
        
        % Draw main figure
        H.FIG=figure('Position',fig_pos, ...
            'Units','Pixel', ...
            'Name','DICOM Loader', ...
            'Numbertitle','off', ...
            'Tag','d_loader_fig', ...
            'Color',GD.col.mainfig, ...
            'Toolbar','none', ...
            'Menubar','none', ...
            'DoubleBuffer','on', ...
            'DockControls','off',...
            'renderer','OpenGL',...
            'CloseRequestFcn',@l_quit,...         
            'Handlevisibility','off');
        
        %                        'KeyPressFcn',@l_KeyPressFcn,...
        
        
        % File Uimenu ---------------------------
        file_h = uimenu('Label','File','Accelerator','F', ...
            'Parent',H.FIG);
        
        H.FILEMENU_OPEN_MRI=uimenu(file_h,'Label','Open &MRI data',...
            'Accelerator','M', ...
            'callback',@l_Open_DATA,...
            'Tag', 'open_MRI',...
            'separator','off');
        
        H.FILEMENU_OPEN_mMRI=uimenu(file_h,'Label','Open M&ultiple MRI datasets',...
            'Accelerator','U', ...
            'callback',@l_Open_DATA,...
            'Tag', 'open_mmri',...
            'separator','off');
        
        H.FILEMENU_OPEN_PET=uimenu(file_h,'Label','Open &PET data',...
            'Accelerator','P', ...
            'callback',@l_Open_DATA,...
            'Tag', 'open_pet',...
            'separator','off');
        
        H.FILEMENU_QUIT = uimenu(file_h,'Label','Exit d_loader','Accelerator','W', ...
            'Callback',@l_quit,'Separator','on');
        
        
        % Main loading panel
        H.LoadPanel=uipanel('parent',H.FIG,...
            'units','pixel',...
            'position',Dat.Pos.LoadPanel,...
            'backgroundcolor',GD.col.frame,...
            'title','Loaders',...
            'visible','on');
        
        H.BTN_MRI=uicontrol('parent',H.LoadPanel,...
            'units','pixel',...
            'position',Dat.Pos.LoadMRI,...
            'style','pushbutton',...
            'string','Load MRI data',...
            'horizontalalignment','center',...
            'callback',@l_Open_DATA,...
            'Tag', 'open_mri',...
            'enable','on',...
            'visible','on');
        
        H.BTN_mMRI=uicontrol('parent',H.LoadPanel,...
            'units','pixel',...
            'position',Dat.Pos.LoadmMRI,...
            'style','pushbutton',...
            'string','Load Multiple MRI datasets',...
            'horizontalalignment','center',...
            'callback',@l_Open_DATA,...
            'Tag', 'open_mmri',...
            'enable','on',...
            'visible','on');
        
        H.BTN_PET=uicontrol('parent',H.LoadPanel,...
            'units','pixel',...
            'position',Dat.Pos.LoadPET,...
            'style','pushbutton',...
            'string','Load PET data',...
            'horizontalalignment','center',...
            'callback',@l_Open_DATA,...
            'Tag', 'open_pet',...
            'enable','on',...
            'visible','on');
        
      
        % Quit btn
        H.BTN_QUIT=uicontrol('parent',H.FIG,...
            'units','pixel',...
            'position',Dat.Pos.Quit,...
            'style','pushbutton',...
            'string','Quit',...
            'horizontalalignment','center',...
            'callback',@l_quit,...
            'enable','on',...
            'visible','on');
        
        
        
        if debug
            debug_h=uimenu('Label','Debug', ...
                'Parent',H.FIG,...
                'enable','on',...
                'callback',@l_debug);
        end
        
    end % l_draw_gui




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure Resize Function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l_ResizeFcn(h,evd)
       %keyboard;
        try
            screen=get(0,'screensize');
            fig_pos=get(H.FIG,'position');
            
            % check if minimum position is requirement is met
            if fig_pos(3)<240 fig_pos(3)=240; end
            if fig_pos(4)<200 fig_pos(4)=200; end
            
            if fig_pos(2)+fig_pos(4)>screen(4)-43
                % top of the figure is outside screen.. bring back.
                fig_pos(2)=screen(4)-fig_pos(4)-43;
            end
            set(H.FIG,'position',fig_pos);
            Dat.Pos.FigPos=fig_pos;
            
            % save new position to preferences..
            Dat.HSiz=fig_pos(3)/screen(3);
            Dat.VSiz=fig_pos(4)/screen(4);
            Dat.HPos=fig_pos(1)/screen(3);
            Dat.VPos=fig_pos(2)/screen(4);
            
            setpref('d_LOADER','HSiz',Dat.HSiz);
            setpref('d_LOADER','VSiz',Dat.VSiz);
            setpref('d_LOADER','HPos',Dat.HPos);
            setpref('d_LOADER','VPos',Dat.VPos);
            
            
            figw=round(Dat.HSiz*screen(3)); % set screen to whatever it was
            figh=round(Dat.VSiz*screen(4)); % set screen to whatever it was
            dx=10;
            dy=10;
            editw=45;
            edith=20;
            textw=120;
            axh=figh/2-2*dy; % axis height;
            axw=figw/2-2*dx; % axis width;
            Dat.Eps=1/axw; % "small nonzero unit"
            loadh=figh-8*dy;
            loadw=figw-2*dx;
            btnh=loadh/3-dy;  % 6 buttons..
            btnw=loadw-1.5*dx;
            
            
            
            % Panels & boxes & others..
            Dat.Pos.LoadPanel=[dx 7*dy loadw loadh];
            Dat.Pos.LoadMRI  =[dx/2 loadh-1*btnh-2*dy btnw btnh];
            Dat.Pos.LoadmMRI =[dx/2 loadh-2*btnh-2*dy btnw btnh];
            Dat.Pos.LoadPET  =[dx/2 loadh-3*btnh-2*dy btnw btnh];
            Dat.Pos.Quit     =[dx*8 dy btnw-10*dx 5*dy];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % resizes
            % Background panel for information and adjusts
            set(H.BTN_QUIT,'position',Dat.Pos.Quit);
            set(H.LoadPanel,'position',Dat.Pos.LoadPanel);
            set(H.BTN_MRI,'position',Dat.Pos.LoadMRI);
            set(H.BTN_mMRI,'position',Dat.Pos.LoadmMRI);
            set(H.BTN_PET,'position',Dat.Pos.LoadPET);
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        catch
            aedes_errordump(lasterror);
        end
    end % l_ResizeFcn






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check preferences
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l_Check_prefs
        
        if ~ispref('d_LOADER','HSiz')
            Dat.HSiz=0.15; % Default window size is 80%
            setpref('d_LOADER','HSiz',Dat.HSiz);
        else
            Dat.HSiz=getpref('d_LOADER','HSiz');
        end
        
        if ~ispref('d_LOADER','VSiz')
            Dat.VSiz=0.2; % Default window size is 80%
            setpref('d_LOADER','VSiz',Dat.VSiz);
        else
            Dat.VSiz=getpref('d_LOADER','VSiz');
        end
        
        if ~ispref('d_LOADER','HPos')
            Dat.HPos=0.1; % Default window size is 80%
            setpref('d_LOADER','HPos',Dat.HPos);
        else
            Dat.HPos=getpref('d_LOADER','HPos');
        end
        
        if ~ispref('d_LOADER','VPos')
            Dat.VPos=0.1; % Default window size is 80%
            setpref('d_LOADER','VPos',Dat.VPos);
        else
            Dat.VPos=getpref('d_LOADER','VPos');
        end
        
        if ~ispref('d_LOADER','LastPath')
            Dat.LastPath=pwd;
            setpref('d_LOADER','LastPath',Dat.LastPath);
        else
            Dat.LastPath=getpref('d_LOADER','LastPath');
            % check that it exists!
            if ~(exist(Dat.LastPath)==7)
                Dat.LastPath=pwd;
            end
        end
        
        
    end % l_Check_prefs


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QUIT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l_quit(h,evd)
        try
            fh=H.FIG;
            
            % Check something should be saved
            %cancel=l_CheckRoiSaved;
            %if cancel
            %  return
            %end
            
            % Save window position
            l_ResizeFcn([],[]);
            % Save preferences
            %l_Check_prefs;
            
            if Dat.To_Aedes
                f=[];
            else
                f=Dat.DATA;
            end
            
            % Clear public variables
            clear Dat H
                       
            % Close GUI
            delete(fh);
            %disp('Bye!');
            
        catch
            aedes_errordump(lasterror);
        end
    end % l_quit




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reset to defaults!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function l_ClearData(h,evd,opt)
        try
            
            % do something perhaps..
            
        catch
            aedes_errordump(lasterror);
        end
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%        LOADERS
%
%
%
%
%
%
%
%



    function l_Open_DATA(h,evd,opt)
        %MRI_LOADER Reads (Siemens) dicom files and parses output for Aedes
        
        % Mikko Nissi 2011, <mikko.nissi@iki.fi>
     
        DEBUG=0;
        
        
        disp(sprintf('Opening: %s',get(h,'tag')));
        
        % Check preferred path
        if exist(Dat.LastPath)==7
            pname=Dat.LastPath;
        else
            pname=pwd;
        end
        
        % prompt for directory to load
        pname=uigetdir(pname,'Select directory containing data');
        if pname==0
            % nothing (cancel) chosen, return
            return;
        end
        
        % quick check..
        if isequal(pname,0)||isempty(pname)
            disp('Something wrong');
            return;
        end
        
        if exist(pname,'dir')
            % all should be ok
        else
            disp('Not a directory');
            return;
        end
       
        % Check if this directory exists in cache
        if isfield(Dat,'DirCache') && ~isempty(Dat.DirCache)
            % something cached
            for ii=1:length(Dat.DirCache)
                if strcmpi(pname,Dat.DirCache(ii).pname)
                    % match found, use this data
                    dirinfo=Dat.DirCache(ii).dirinfo;
                    % Save path into prefs
                    Dat.LastPath=pname;
                    setpref('d_LOADER','LastPath',Dat.LastPath);
                    break;
                end
            end
        end
        if ~exist('dirinfo','var')
          
            if(DEBUG)
                disp('dirinfo not exist');
            end
          
          % Ok - load stuff manually!
        
          % Apparently ok, save path into prefs
          Dat.LastPath=pname;
          setpref('d_LOADER','LastPath',Dat.LastPath);
        
        
          % dir all files..
          flist3=dir(pname);
        
          flist=flist3;
          % and remove "." and ".."
          flist=flist(3:end);
          flist3=flist3(3:end);
          %keyboard;
        
          % 1. Remove directories and leave files!
          ii=1;
          while ii<=length(flist)
            if flist(ii).isdir
              flist(ii)=[];
            elseif strfind(flist(ii).name,Dat.d_indexfile)   % don't load d_loader's own index..
              flist(ii)=[];
            elseif flist(ii).name(1)=='.'  % if just a dot-file, these are not hidden on mac/unix. this rule might work for all..
              flist(ii)=[];
            else
              ii=ii+1;
            end
          end
          
          % 2. Remove files and leave directories!
          ii=1;
          while ii<=length(flist3)
            if ~flist3(ii).isdir
              flist3(ii)=[];
            else
              ii=ii+1;
            end
          end
          
          if isempty(flist) && isempty(flist3)
            % apparently nothing found..??
            disp('No files found, check your directory or modify the loader program!');
            return;
          end
          
        
          %%%% check subdirs if necessary
          if ~isempty(flist3)>0
            tmplist=dir(fullfile(pname,flist3(1).name));
            for kk=1:length(tmplist)
              tmplist(kk).name=[flist3(1).name filesep tmplist(kk).name];
            end
            flist2=tmplist;
            if length(flist3)>1
              for ii=2:length(flist3)
                tmplist=dir(fullfile(pname,flist3(ii).name));
                for kk=1:length(tmplist)
                  tmplist(kk).name=[flist3(ii).name filesep tmplist(kk).name];
                end
                flist2(end+1:end+length(tmplist))=tmplist;
              end
            end
            
            % Clean non-dicom, indexfiles, '.' and '..' from filelist:
            ii=1;
           hwb=waitbar(0,'Please wait, checking files..');
           modder=floor(length(flist2)/100);
           
           
            while ii<=length(flist2)
              if strfind(flist2(ii).name,Dat.d_indexfile)
                flist2(ii)=[];
              elseif strfind(flist2(ii).name,[filesep '.'])
                flist2(ii)=[];
              elseif strfind(flist2(ii).name,[filesep '..'])
                flist2(ii)=[];
              elseif ~isdicom(fullfile(pname,flist2(ii).name))
                  flist2(ii)=[];
              else
                ii=ii+1;
              end
              if mod(ii,modder)==0
                  waitbar(ii/length(flist2),hwb);
              end
            end
            delete(hwb);
            
            
            % cleaned! take into use!
            flist(end+1:end+length(flist2))=flist2;
          end
          flist_loaded=flist;
          temp_flist=flist;
        
          %keyboard;
        
          % Read dicom info from all files, OR, read index-matrix
          clear dicomheads;
          % reset number of removed files..
          removed_files=[];
          if (exist(fullfile(pname,Dat.d_indexfile))==2)
            hwb=waitbar(0,'Please wait, loading previously saved information');
            % load index
            load(fullfile(pname,Dat.d_indexfile));
            waitbar(1,hwb,'Please wait, loading previously saved information');
            delete(hwb);
          end
        
          if(DEBUG) disp('comparing dir-list to loaded file');end

          
          % Check if files need to be removed from the loaded flist
          if ~isempty(removed_files)
            temp_flist(removed_files)=[];
          end
          
       %   keyboard;
        
          % Check that index is ok or exists: if no, load data
          if ~exist('dicomheads','var')||length(dicomheads)~=length(temp_flist)
              if(DEBUG) disp(sprintf('dicomheads exist: %d',exist('dicomheads','var'))); end
              if(DEBUG) if (exist('dicomheads','var')) disp(sprintf('length dicomheads: %d',length(dicomheads)));end;end
              if(DEBUG) disp(sprintf('length flist_load: %d',length(temp_flist)));end
              %if(DEBUG) keyboard; end;
              
            clear dicomheads; % Just in case it did exist but size was changed..!             
            flist=flist_loaded;
            modder=floor(length(flist)/100);
            hwb=waitbar(0,sprintf('Please wait, checking files: 1/%1.0d',length(flist)));
            ii=1;rci=0;
            removed_files=[];
            saveheads={};
            %saveheads{length(flist)}=[];
            while ii<=length(flist)              
              if mod(ii,modder)==0
                waitbar(ii/length(flist),hwb,sprintf('Please wait, checking files: %1.0d/%1.0d',ii,length(flist)));
              end
              try
                rci=rci+1;
                dicomheads{ii}=dicominfo(fullfile(pname,flist(ii).name));
                ii=ii+1;
              catch
                % did not succeed, remove this file from list and move on
                flist(ii)=[];
                % save the number of files to be removed from the list.. These will
                % be removed later on again if reloaded.. Note that this is likely
                % to FAIL. Ie. if n files are removed and n files added, this is
                % almost guaranteed to blow up..
                removed_files=[removed_files rci];
              end
              % Keep only something..
              KeepFields={'SeriesNumber','SeriesType','SequenceName','ProtocolName','ImageType','FlipAngle','RescaleIntercept','RescaleSlope','SliceLocation','InstanceNumber','EchoTime','AcquisitionNumber'};
              if exist('dicomheads','var') && ~isempty(dicomheads{ii-1})
                  for kk=1:length(KeepFields)
                      %                keyboard;
                      if isfield(dicomheads{ii-1},KeepFields{kk})
                          eval(sprintf('saveheads{ii-1}.%s=dicomheads{ii-1}.%s;',KeepFields{kk},KeepFields{kk}));
                      end
                  end
              end
            end
            % switch to keeper info:
            dicomheads=saveheads;
            delete(hwb);
          end


            %        keyboard;
          % Sort file by SeriesNumber field
          for ii=1:length(dicomheads)
              % Sorting by SeriesNumber, which should be used normally:
              tmp=dicomheads{ii}.SeriesNumber;
              
              % Sorting by FILENAME:
              % should not normally be used.
              % tmp=ii;
              
            if ~isempty(tmp)
              if isnumeric(tmp)
                sortnum(ii)=tmp;
              elseif ischar(tmp)
                sortnum(ii)=str2num(tmp);
              else
                sortnum(ii)=0;
              end
            else
              sortnum(ii)=0;
            end
          end
          
          %%%%%%%%%
          % UNCOMMENT THESE!!!
          
          [sortval,sortndx]=sort(sortnum);
          dicomheads=dicomheads(sortndx);
          flist=flist(sortndx);
        
          %
          %
          %%%%%%
          
          %%% CHECK LISTING AND IF SEVERAL SEQUENCES FOUND, give choice!
          %dprops={}; %empty filename-property-holder
          %for nn=1:length(flist)
          %  dprops(nn,:)=regexp(flist(nn).name,'[^\.]*','match');
          %end
          %for nn2=1:size(dprops,2)
          %  tmp=str2num(dprops{1,nn2});
          %  if ~isempty(tmp)
          %    break;
          %  end
          %end
        
          % vectorize first numeric part of the file name
          %ltmp=cell2mat(dprops(:,nn2));
          %listvec=str2num(ltmp);
        
          % ASSUME THIS INDICATES "NEW" SERIE in the series..
          fstarts=[1 find(diff(sortval)~=0)+1];
          fstarts=[fstarts length(flist)+1]; % add number of the last measurement..
          
          %keyboard;
          % Store stuff in cache
          dirinfo.fstarts=fstarts;
          dirinfo.dicomheads=dicomheads;
          dirinfo.flist=flist;
          dirinfo.removed_files=removed_files;
          
          % Store to glob var
          if(isfield(Dat,'DirCache'))
            Dat.DirCache(end+1).dirinfo=dirinfo;
            Dat.DirCache(end).pname=pname;
          else
            Dat.DirCache.dirinfo=dirinfo;
            Dat.DirCache.pname=pname;
          end
            
          % save data
          hwb=waitbar(1,'Saving information for later');
          save(fullfile(pname,Dat.d_indexfile),'dicomheads','removed_files','flist');
          delete(hwb);

        end
 
        %keyboard;
        % Now all info either taken from cache or reloaded, stuff back to
        % vars
        fstarts=dirinfo.fstarts;
        dicomheads=dirinfo.dicomheads;
        flist=dirinfo.flist;
        removed_files=dirinfo.removed_files;
        
        %keyboard;
        
        % Sort data further depending on what's being loaded!
        switch get(h,'tag')
            case 'open_pet'
                clear istr istr2;
                for jj=1:length(fstarts)-1
                    itmp=dicomheads{fstarts(jj)};
                    istr2{jj}=sprintf('%03.0f: files %4.0f - %4.0f (n = %3.0f), series: %15s',jj,fstarts(jj),fstarts(jj+1)-1,fstarts(jj+1)-fstarts(jj),itmp.SeriesType);
                    istr{jj}=sprintf('%03.0f: %3.0f files, seq: %15s, series: %25s ',jj,fstarts(jj+1)-fstarts(jj),itmp.SeriesType);
                    %  disp(istr2{jj});
                end
                
          otherwise
                clear istr istr2;
                for jj=1:length(fstarts)-1
                    itmp=dicomheads{fstarts(jj)};
                    if isfield(itmp,'ProtocolName')
                        if isfield(itmp,'SequenceName') seqNam=itmp.SequenceName; else seqNam='';end
                        if isfield(itmp,'ProtocolName') protNam=itmp.ProtocolName; else protNam='';end
                        if isfield(itmp,'SeriesNumber') serNum=itmp.SeriesNumber; else serNum='';end
                        if isfield(itmp,'ImageType') imType=itmp.ImageType; else imType='';end
                        istr2{jj}=sprintf('%8s/ %03.0f: files %4.0f - %4.0f (n = %3.0f), seq: %15s, protocol: %25s, SeriesNumber: %1.0f ',fileparts(flist(fstarts(jj)).name),jj,fstarts(jj),fstarts(jj+1)-1,fstarts(jj+1)-fstarts(jj),seqNam,protNam,serNum);
                        istr{jj}=sprintf('%8s/ %03.0f: %3.0f files, seq: %15s, protocol: %25s, imagetype: %25s ',fileparts(flist(fstarts(jj)).name),jj,fstarts(jj+1)-fstarts(jj),seqNam,protNam,imType);
                    else
                        if isfield(itmp,'SequenceName') seqNam=itmp.SequenceName; else seqNam='';end
                        if isfield(itmp,'ProtocolName') protNam=itmp.ProtocolName; else protNam='';end
                        if isfield(itmp,'SeriesNumber') serNum=itmp.SeriesNumber; else serNum='';end
                        if isfield(itmp,'ImageType') imType=itmp.ImageType; else imType='';end
                        istr2{jj}=sprintf('%8s/ %03.0f: files %4.0f - %4.0f (n = %3.0f), seq: %15s, SeriesNumber: %1.0f ',fileparts(flist(fstarts(jj)).name),jj,fstarts(jj),fstarts(jj+1)-1,fstarts(jj+1)-fstarts(jj),seqNam,serNum);
                        istr{jj}=sprintf('%8s/ %03.0f: %3.0f files, seq: %15s, imagetype: %25s ',fileparts(flist(fstarts(jj)).name),jj,fstarts(jj+1)-fstarts(jj),seqNam,imType);
                    end
                    %  disp(istr2{jj});
                end
        end
        
        
        % If loading multiple data sets, give multiple selection!
        if strcmp(get(h,'tag'),'open_mmri')
            % Give selection box FOR MULTIPLE FILES!!
            [s,v] = d_filedlg('PromptString',Dat.title,'ListString',istr,'listsize',[700,300],'selectionmode','multiple');
        else
            % Give selection box FOR SINGLE FILE:
            [s,v] = d_filedlg('PromptString',Dat.title,'ListString',istr,'listsize',[700,300],'selectionmode','single');
        end
        
   %     keyboard;
        
        % adjust filelist accordingly!
        if isempty(s) | v==0
            % nothing selected
            return;
            
        end
        
        % Set a no-flag for multiple 4-D / 5-D loading...
        fiveD_flag=0;
        
        % loop selection and load files into same var
        for ii=1:length(s)
            % remove all extras:
            filelist=flist(fstarts(s(ii)):fstarts(s(ii)+1)-1);
            headlist=dicomheads(fstarts(s(ii)):fstarts(s(ii)+1)-1);
           
            % Load files:
            [tmpdat,dcm_headertmp]=l_fileprocess(pname,filelist,headlist,h);
            if ii==1
                infotmp(ii,:)=dcm_headertmp;
            else
                % Second set of files -- attempt to merge next set of DICOM headers
                infotmp=l_Merge_Dicom_Info(infotmp,dcm_headertmp);
            end

            %keyboard;
            
            if length(s)>1               

                if ndims(tmpdat)==4
                    % there are too many dimensions!!
                    % Attempt loading all and re-shuffling data to extended
                    % 4-D matrix. If that also fails, then complain
                    if fiveD_flag==0
                        % warn on the first file only..
                        disp('Warning - loading of multiple 4D datasets attempted.');
                        disp(' \__Attempting to load and concatenate.');                        
                        disp('      \__Fails if volumes are not of the same size!');
                    end
                    % reset 4-D / 5-D flag
                    fiveD_flag=1;
                    tmp2(:,:,:,:,ii)=tmpdat;
                else
                    % Ok, loading multiple 2D or 3D sets.. (actually may
                    % not work with multiple 2Ds anyways)
                    tmp2(:,:,:,ii)=tmpdat;
                end
            else
                % single file, just load whatever is there..
                tmp2=tmpdat;
            end
        end
        
        %keyboard;
        
        % Check if multiple 4-D sets were loaded and attempt to reshape
        if fiveD_flag
            tmp2=reshape(tmp2,[size(tmp2,1) size(tmp2,2) size(tmp2,3) 2*size(tmp2,4)]);
        end
        
        % sort files by flip angle if multiple (DISABLED FOR MULTIPLE 4-D SETS --- MAY
        % NOT WORK..
        if ~fiveD_flag && strcmp(get(h,'tag'),'open_mmri') && length(s)>1
            for kk=1:length(s)
                sindex(kk)=dicomheads{fstarts(s(kk))}.FlipAngle;
            end
            %  keyboard;
            [val,ndx]=sort(sindex);
            tmp2=tmp2(:,:,:,ndx);
        else
            % nothing here..
        end
        
        
        %keyboard;
        % load aedes with this data!
        DATA.FTDATA=tmp2;
        DATA.HDR.dicominfo=infotmp;
        DATA.HDR.FileHeader=infotmp;
        DATA.HDR.fpath = [pname filesep];
        DATA.HDR.fname = filelist(1).name; % put the first name in
        DATA.HDR.FileFormatName = 'DICOM';
        DATA.DataFormat = 'dcm';
        
        % Structure in order, open Aedes:
        if Dat.To_Aedes
            aedes(DATA);
        else
            % No aedes, but output requested.. Stuff data and quit!
            Dat.DATA=DATA;
            l_quit(h,evd);
        end
        

        
    end




%%
%%%%%%%%%%%%%%%%%%
% SUBFUNCTION FOR LOADING FILES


    function f=l_Merge_Dicom_Info(dcm_1,dcm_2)
        % attempt to concatenate dcm_2 after dcm_1, cross-matching whatever
        % fields necessary to make the headers compatible
        
        fields_1=fieldnames(dcm_1(1)); % just take the first since the group is assumed same anyways
        fields_2=fieldnames(dcm_2(1)); % just take the first since the group is assumed same anyways
        
        % Go through dcm_1 and add any missing fields to dcm_2
        %        rrtmp=infotmp(1);
        %keyboard;
        
        for xxk=1:length(fields_1)
            if ~isfield(dcm_2,fields_1{xxk})
                %disp(fields_1{xxk});
                % set the same field to dcm_2, taking value from dcm_1:
                for ii=1:length(dcm_2)
                    % It is assumed there is the same number of images in
                    % both sets. Otherwise something else will fail...
                    dcm_2(ii).(fields_1{xxk})=dcm_1(ii).(fields_1{xxk});
                end
            end
        end
        
        % Repeat the process the other way
        for xxk=1:length(fields_2)
            if ~isfield(dcm_1,fields_2{xxk})
                %disp(fields_2{xxk});
                % set the same field to dcm_1, taking value from dcm_2:
                for ii=1:length(dcm_2)
                    % It is assumed there is the same number of images in
                    dcm_1(ii).(fields_2{xxk})=dcm_2(ii).(fields_2{xxk});
                end
            end
        end

        % Concatenate dcms - requires re-ordering the fields to make the
        % structures concate
        f=orderfields(dcm_1);
        f(2,:)=orderfields(dcm_2);
        return;
        
    end

% End l_Merge_Dicom_Info(dcm_1,dcm_2)



    function [f,infotmp]=l_fileprocess(pname,flist,dicomheads,h)
        
        
        % check image size from first image:
        hwb=waitbar(0,sprintf('Loading file: 1/%1.0d',length(flist)));
        ctmp=dicomread(fullfile(pname,flist(1).name));
        %        itmp=dicomheads{1};
        % read data from file - this info is not available anymore..
        itmp=dicominfo(fullfile(pname,flist(1).name));
        tmp=double(zeros(size(ctmp,1),size(ctmp,2),length(flist),class(ctmp)));
        infotmp=itmp;
        infotmp(length(flist))=itmp;
        
        % Check if scaled:
        if isfield(infotmp(1),'RescaleIntercept') && isfield(infotmp(1),'RescaleSlope')
            ctmp=infotmp(1).RescaleIntercept+double(ctmp).*infotmp(1).RescaleSlope;
            disp(sprintf('Scaling image data: y = %1.3f * x + %1.3f',infotmp(1).RescaleSlope,infotmp(1).RescaleIntercept));
        else
            % rescale not found, SET IT!! To help process files properly
            % later..
            infotmp(1).RescaleIntercept=0;
            infotmp(1).RescaleSlope=1;            
        end
        
        %keyboard;       
        tmp(:,:,1)=ctmp;
        
        count=1;
        ndx=1;
        % start loading!
        for ii=2:length(flist)
            waitbar(ii/length(flist),hwb,sprintf('Loading file: %1.0d/%1.0d',ii,length(flist)));
            ctmp=dicomread(fullfile(pname,flist(ii).name));
            %            itmp=dicomheads{ii}; %dicominfo(fullfile(pname,flist(ii).name));
            % Again, read header info from file..
            itmp=dicominfo(fullfile(pname,flist(ii).name));
            if isequal(size(tmp(:,:,1)),size(ctmp))
                %   keyboard;
                count=count+1;
                
                % Check if scaled:
                if isfield(itmp,'RescaleIntercept') && isfield(itmp,'RescaleSlope')
                    ctmp=itmp.RescaleIntercept+double(ctmp).*itmp.RescaleSlope;
              %      disp(sprintf('Scaled data = %1.3f * x + %1.3f',itmp.RescaleSlope,itmp.RescaleIntercept));
              %      keyboard;
                else
                    % rescale not found, SET IT!! To help process files properly
                    % later..
                    itmp(1).RescaleIntercept=0;
                    itmp(1).RescaleSlope=1;
                end
                
                tmp(:,:,count)=ctmp;
               
                % Apparently dicom header fields may change within one sequence, do
                % something to avoid confusion: take ONLY those fields present in the
                % first image and dump whatever is added or missing later on..
                fields=fieldnames(infotmp(1));
                rrtmp=infotmp(1);
                for xxk=1:length(fields)
                    if isfield(itmp,fields{xxk})
                        rrtmp=setfield(rrtmp,fields{xxk},getfield(itmp,fields{xxk}));
                    %else
                    %    rrtmp=setfield(rrtmp,fields{xxk},getfield(itmp,fields{xxk}));
                    end
                    infotmp(count)=rrtmp;
                end
                %    keyboard;
                ndx(count)=ii;
            else
                % size does not match!
                %    keyboard;
                disp(sprintf('Size of file %s does not match the first of the file, skipping!',flist(ii).name));
                %infotmp(count+1)=[];
            end
        end
        delete(hwb);
        
        %keyboard;
        
        % loaded, determine number of slices
        % scrap excess stuff:
        tmp=tmp(:,:,ndx);
        infotmp=infotmp(ndx);
        ntmp1={};
        ntmp1b={};
        if isfield(infotmp(1),'SliceLocation')
            [ntmp1{1:count,1}]=deal(infotmp.SliceLocation);
            ntmp2=cell2mat(ntmp1); % put slice locations into single var
        end
        if isfield(infotmp(1),'InstanceNumber')           
            [ntmp1b{1:count,1}]=deal(infotmp.InstanceNumber);
            ntmp2b=cell2mat(ntmp1b); % put instance numbers into single var
        else
            ntmp2b=ntmp2;  % in case instance numbers not found, use slice for both
        end
        nslices=length(find(diff(sort(ntmp2))))+1;
        
        
        
      %  keyboard;
        
        
        % Sort data further depending on what's being loaded!
        %switch get(h,'tag')
        %    case {'open_mri','open_mmri'}
        %        % Try sorting the data few times:
        %        %  1. by slice location
        %        [ssloc,sortndx]=sort(ntmp2);
        %        tmp=tmp(:,:,sortndx);
        %        infotmp=infotmp(sortndx);
        %
        %        %  2. by echo times!
        %        echotmp1={};
        %        [echotmp1{1:count,1}]=deal(infotmp.EchoTime);
        %        echotmp2=cell2mat(echotmp1); % put slice locations into single var
        %        %   sort with echo times
        %        [ssloc,sortndx2]=sort(echotmp2);
        %        tmp=tmp(:,:,sortndx2);
        %        infotmp=infotmp(sortndx2);
        %
        %        %  3. by acquisitionnumber
        %        acqnumtmp={};
        %        [acqnumtmp{1:count,1}]=deal(infotmp.AcquisitionNumber);
        %        acqnumtmp2=cell2mat(acqnumtmp);
        %        % sort
        %        [ssloc,sortndx2]=sort(acqnumtmp2);
        %        tmp=tmp(:,:,sortndx2);
        %        infotmp=infotmp(sortndx2);
        %
        %    case {'open_pet'}
        %        % Sort data along slice locations!
        %        [ssloc,sortndx]=sort(ntmp2);
        %        tmp=tmp(:,:,sortndx);
        %        infotmp=infotmp(sortndx);
        %
        %    otherwise
        
                       
        % Sort data along instance numbers first! (same as slice if this
        % field doesn't exist)
    %    [ssloc,sortndx]=sort(ntmp2b);
    %    tmp=tmp(:,:,sortndx);
    %    infotmp=infotmp(sortndx);
                
        
        % Sort data along slice locations!
        [ssloc,sortndx]=sort(ntmp2);
        tmp=tmp(:,:,sortndx);
        infotmp=infotmp(sortndx);
        
        % TRY SORTING MORE!
        %  2. by echo times!
        if(isfield(infotmp,'EchoTime'))
            echotmp1={};
            [echotmp1{1:count,1}]=deal(infotmp.EchoTime);
            echotmp2=cell2mat(echotmp1); % put slice locations into single var
            %   sort with echo times
            if ~isempty(echotmp2)
                [ssloc,sortndx2]=sort(echotmp2);
                tmp=tmp(:,:,sortndx2);
                infotmp=infotmp(sortndx2);
            end
        end
        
        %  3. by acquisitionnumber
        if(isfield(infotmp,'AcquisitionNumber'))
            acqnumtmp={};
            [acqnumtmp{1:count,1}]=deal(infotmp.AcquisitionNumber);
            acqnumtmp2=cell2mat(acqnumtmp);
            % sort
            if ~isempty(acqnumtmp2)
                [ssloc,sortndx2]=sort(acqnumtmp2);
                tmp=tmp(:,:,sortndx2);
                infotmp=infotmp(sortndx2);
            end
        end
        %end
        
        %  4. by acquisitiontime
        if(isfield(infotmp,'AcquisitionTime'))
            acqnumtmp={};
            [acqnumtmp{1:count,1}]=deal(infotmp.AcquisitionTime);
            acqnumtmp2=str2num(cell2mat(acqnumtmp));
            % sort
            if ~isempty(acqnumtmp2)
                [ssloc,sortndx2]=sort(acqnumtmp2);
                tmp=tmp(:,:,sortndx2);
                infotmp=infotmp(sortndx2);
            end
        end
        %end

        
        
        
        % reshape data according to this;
        f=reshape(tmp,[size(tmp,1),size(tmp,2),nslices,count/nslices]);
        
        % remove if singleton dimensions added.... 
        f=squeeze(f);
    end
















%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

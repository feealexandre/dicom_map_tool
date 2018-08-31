function varargout = dicom_map_tool(varargin)
% DICOM_MAP_TOOL MATLAB code for dicom_map_tool.fig
%      DICOM_MAP_TOOL, by itself, creates a new DICOM_MAP_TOOL or raises the existing
%      singleton*.
%
%      H = DICOM_MAP_TOOL returns the handle to a new DICOM_MAP_TOOL or the handle to
%      the existing singleton*.
%
%      DICOM_MAP_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICOM_MAP_TOOL.M with the given input arguments.
%
%      DICOM_MAP_TOOL('Property','Value',...) creates a new DICOM_MAP_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dicom_map_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dicom_map_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Revisions

%   02/2018, Felipe Alexandre
%      - Added button to load multiples niftis files of T1r and T2r. The 
%        recognition is made by the name of the files, this way you can load
%        differents kinds of files names, in case that it not contain T1r or T2r 
%        inside the name, will be skip. Also compares the time sequence setted 
%        with time sequence of the files. The output is saving automatically with
%        '_map' in the end of name.
%      - Added raff input options to choose between load 1 file with (+Z) and
%        (-Z) or load 2 files (+Z and -Z).
%      - Added button to convert 2 raff files (+Z and -Z) in 1 raff files with 
%        both. The output is saving automatically using the same name of the inv
%        files, but with '_plu' in the end of the name.
%      - Added button to load multiples niftis files of raffs. The user can choose
%        the type of inputs, before load multiples files. In case of choose 1 input
%        files, the script will calculate for file with "raff" or "RAFF" in the name
%        In case of choose 2 input files, the name of the files need to follow the rules. 
%        Example:
%        name of a +Z file : "2219A_15s_gre_prep_raff3_plu_tr5_tra.nii"
%        name of a -Z file : "2219A_16s_gre_prep_raff3_inv_tr5_tra.nii"
%        The software will delete since the first '_' until the second '_'
%        because is for numeration of files, and delete the '_plu' and '_inv'
%        after it, to join the files and calculate it must has the same name. 
% 05/2016, Mikko Nissi
%    - Added 3-D volume mask-generation to 4-D data instead of slice-based
%      which obviously may fail with noise-only slices..
%    - Added range limits to both T1r & MT fitting. The maximum/minimum of
%      saved relaxation time values are now restricted to [0 3000].
%
% 04/2016, Mikko Nissi
%    - Added Nifti-saving options
%    - Nifti-save will check if the input was Nifti and in that case uses
%      the original header to save any coordinate and other information
%    - MAY break scaling of data values if Nifti-header has such information
%    - Not thoroughly tested yet.
%
% 03/2016, Mikko Nissi
%    - Fixed a bunch of bugs related to handling single-slice DICOM data
%    - Added recognition of single-slice data in mask generation that works for RAFF
%      data (may break with single-slice T1rho though..)
%
% 03/2016, Mikko Nissi
%    - Now changes working directory to the loaded first RAFF/MT dataset.
%      This should help loading the next dataset correctly, but may cause
%      other confusion.
%
% 02/2016, Mikko Nissi
%    - Added Nifti input / reading functionality 
%    - Added manual input for tSL/TE/etc
%    - Corrected couple of bugs here and there..
%
% 01/2014, Mikko Nissi
%    - Fixed a critical bug forcing use of MT timings for RAFF fitting.
%
% (c) Mikko Nissi 2013-2014



% Edit the above text to modify the response to help dicom_map_tool

% Last Modified by GUIDE v2.5 02-Feb-2018 15:48:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dicom_map_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @dicom_map_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dicom_map_tool is made visible.
function dicom_map_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dicom_map_tool (see VARARGIN)

% Choose default command line output for dicom_map_tool
handles.output = hObject;
handles.DEBUG=0;
handles.NiftiInput=0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dicom_map_tool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dicom_map_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_quit.
function btn_quit_Callback(hObject, eventdata, handles)
% hObject    handle to btn_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);





% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2






%%%%%%%%%%%%%%%%%%%%%  LOADING BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%
%

% --- Executes on button press in btn_loadNifti_T1_for_RAFFMT.
function btn_loadNifti_T1_for_RAFFMT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_T1_for_RAFFMT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('T1 not implemented, aborting...');
return;



% --- Executes on button press in btn_loadNifti_T1_for_T12r.
function btn_loadNifti_T1_for_T12r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_T1_for_T12r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('T1 not implemented, aborting...');
return;


% --- Executes on button press in btn_load_t1_for_mt.
function btn_load_t1_for_mt_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_t1_for_mt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T1 map data to handles structure, run d_loader with SINGLE file
% selection:

disp('T1 not implemented, aborting...');
return;

handles.T1MAP=d_loader('openmode',1,'title','Please choose data for T1 relaxation time map:');

% ----------TODO----------CHECK IF T1MAP ok..? 

if ~isempty(handles.T1MAP)
    % Enable 3-parameter-fitting with T1
    set(handles.radiobutton13,'enable','on');
    set(handles.radiobutton14,'enable','on');
end


% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_load_mt.
function btn_load_mt_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_mt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load MT data to handles structure, run d_loader with MULTIPLE file
% selection:
handles.SOURCE=d_loader('openmode',2,'title','Please choose +Z and -Z source data for MT:');

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_mt,'enable','on');
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=2;
handles.NiftiInput=0;
handles.SOURCETYPE='mt_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_loadNifti_MT.
function btn_loadNifti_MT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_MT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load RAFF data to handles structure, run d_loader with MULTIPLE file
% selection:


% Load MT data to handles structure FROM NIFTI SOURCE WITH TWO READS
uiwait(msgbox('Please open the first MT (+Z) dataset from the next dialog!','Note!','modal'));
tmp1=aedes_read_nifti;
    
% CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
% NOTE // THIS MAY CAUSE OTHER CONFUSION...
if ~isempty(tmp1) && isfield(tmp1,'HDR') && isfield(tmp1.HDR,'fpath') && ~isempty(tmp1.HDR.fpath)
    cd(tmp1.HDR.fpath);
end

uiwait(msgbox('Now please open the second MT (-Z) dataset (or none) from the next dialog!','Note!','modal'));
tmp2=aedes_read_nifti;

% Append the second dataset to the first:
if ~isempty(tmp2)
    tmp1.FTDATA=cat(ndims(tmp1.FTDATA),tmp1.FTDATA,tmp2.FTDATA);
end

handles.SOURCE=tmp1;


if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_mt,'enable','on');
    handles.SOURCE.FTDATA=double(handles.SOURCE.FTDATA);
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=2;
handles.NiftiInput=1;
handles.SOURCETYPE='mt_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);






% --- Executes on button press in btn_load_raff.
function btn_load_raff_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_raff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load RAFF data to handles structure, run d_loader with MULTIPLE file
% selection:
handles.SOURCE=d_loader('openmode',2,'title','Please choose +Z and -Z source data for RAFF:');

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_mt,'enable','on');
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=2;
handles.NiftiInput=0;
handles.SOURCETYPE='raff_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_loadNifti_RAFF.
function btn_loadNifti_RAFF_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_RAFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load RAFF data to handles structure FROM NIFTI SOURCE WITH TWO READS
uiwait(msgbox('Please open the first RAFF (+Z) dataset from the next dialog!','Note!','modal'));
tmp1=aedes_read_nifti;

% CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
% NOTE // THIS MAY CAUSE OTHER CONFUSION...
if ~isempty(tmp1) && isfield(tmp1,'HDR') && isfield(tmp1.HDR,'fpath') && ~isempty(tmp1.HDR.fpath)
    cd(tmp1.HDR.fpath);
end

uiwait(msgbox('Now please open the second RAFF (-Z) dataset (or none) from the next dialog!','Note!','modal'));
tmp2=aedes_read_nifti;

% Append the second dataset to the first:
if ~isempty(tmp2)
    tmp1.FTDATA=cat(ndims(tmp1.FTDATA),tmp1.FTDATA,tmp2.FTDATA);
end

handles.SOURCE=tmp1;


if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_mt,'enable','on');
    handles.SOURCE.FTDATA=double(handles.SOURCE.FTDATA);
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=2;
handles.NiftiInput=1;
handles.SOURCETYPE='raff_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_load_t1r.
function btn_load_t1r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_t1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T1rho data to handles structure, run d_loader with SINGLE file
% selection:
handles.SOURCE=d_loader('openmode',1,'title','Please choose source data for T1r:');

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_t1r,'enable','on');
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=1;
handles.NiftiInput=0;
handles.SOURCETYPE='t1rho_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);


% --- Executes on button press in btn_loadNifti_T1r.
function btn_loadNifti_T1r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_T1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T1rho data to handles structure FROM NIFTI SOURCE
handles.SOURCE=aedes_read_nifti;

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_t1r,'enable','on');
    handles.SOURCE.FTDATA=double(handles.SOURCE.FTDATA);
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=1;
handles.NiftiInput=1;
handles.SOURCETYPE='t1rho_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);




% --- Executes on button press in btn_load_t2r.
function btn_load_t2r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_t2r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T2rho data to handles structure, run d_loader with SINGLE file
% selection:
handles.SOURCE=d_loader('openmode',1,'title','Please choose source data for T2r:');

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_t1r,'enable','on');
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=1;
handles.NiftiInput=0;
handles.SOURCETYPE='t2rho_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);


% --- Executes on button press in btn_loadNifti_T2r.
function btn_loadNifti_T2r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadNifti_T2r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T2rho data to handles structure FROM NIFTI SOURCE
handles.SOURCE=aedes_read_nifti;

if ~isempty(handles.SOURCE)
    % Enable fitting-button
    set(handles.btn_fit_t1r,'enable','on');
    handles.SOURCE.FTDATA=double(handles.SOURCE.FTDATA);
end

% set a tag for identifying MT/RAFF vs. T12rho
handles.LR=1;
handles.NiftiInput=1;
handles.SOURCETYPE='t2rho_map_';

% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_load_t1_for_t1r.
function btn_load_t1_for_t1r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_load_t1_for_t1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load T1 map data to handles structure, run d_loader with SINGLE file
% selection:
handles.T1MAP=d_loader('openmode',1,'title','Please choose data for T1 relaxation time map:');

% ----------TODO----------CHECK IF T1MAP ok..? 

if ~isempty(handles.T1MAP)
    % Enable 3-parameter-fitting with T1
    set(handles.radiobutton4,'enable','on');
end

% Stuff back to GUI for later use..
guidata(hObject,handles);







% --- Executes on button press in check_savenifti.
function check_savenifti_Callback(hObject, eventdata, handles)
% hObject    handle to check_savenifti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_savenifti

% Don't do anything right now.. Set on per default to always save a Nifti..
% --- Executes on button press in btn_loadultiNifti_T1r.



function btn_loadMultiNifti_T1r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadMultiNifti_T1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load T1rho data to handles structure FROM NIFTI SOURCE

%Load multiples files, f_name is variable with all the names of
%the files,f_path with paths.
%aedes_juigetfiles is a function from aedes.
[f_name,f_path,f_index] = aedes_juigetfiles;


n_files = size(f_name,2); %N_files is the number of files loaded

% Go through all the file selecting one by one to save in handle structure
    for n_files = 1 : n_files
        % Check if the user loaded same data, prevent debugs
        if(f_index ~= 0)
        %Set datas on handles structure, making a vector of datas. 
        handles.SOURCE(n_files) = aedes_data_read([f_path{n_files},f_name{n_files}]);

            if ~isempty(handles.SOURCE)
                % Enable fitting-button
                set(handles.btn_fit_t1r,'enable','on');
                handles.SOURCE(n_files).FTDATA=double(handles.SOURCE(n_files).FTDATA);
                % Disable the output-options
                set(handles.radio_saveopen_DICOM1,'enable','off'); 
                set(handles.radio_saveopen_AEDES1,'enable','off');
            end
        end
    end


% Stuff back to GUI for later use..
guidata(hObject,handles);



% --- Executes on button press in btn_loadMultiNifti_raff.
function btn_loadMultiNifti_raff_Callback(hObject, eventdata, handles)
% hObject    handle to btn_loadMultiNifti_raff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Getting the type of inputs
input_mask = [get(handles.radiobutton15,'Value'),...
            get(handles.radiobutton16,'Value')];
i = 1;
p = 1;
c = 1;
control3 = 1;

switch find(input_mask)
    % Case 1 the input files will be loaded RAFF (+Z) and (-Z) in one nifti file 
    case 1
        uiwait(msgbox('Please open multiples files with RAFF (+Z) and RAFF (-Z) in 1 nifti file from the next dialog!','Note','modal'));
        % Load RAFF data to handles structure FROM NIFTI SOURCE WITH
        % MULTIPLES READS OF ONE FILE WITH (+Z) e (-Z).
        [f_name,f_path,f_index] = aedes_juigetfiles;
        n_files = size(f_name,2);
        control1= n_files;
        % Check if the user loaded same data, prevent debugs
        if(f_index ~= 0)
        % Go through all the file saving one by one in handles structure
        for control1 = 1:control1
                    % Skiping the files that has '_map' on name (output files).
                    if(contains(f_name{control1},'_map')==1)
                        disp(['The file ' f_name{control1} ' was skiped']);
                    continue;
                    end
                    % Checking if the files has 'raff'/'RAFF' on name.
                    if contains(f_name{control1},'raff')==1 || contains(f_name{control1},'RAFF')==1
                        handles.SOURCE(control1) = aedes_data_read([f_path{control1},f_name{control1}]);
                    else
                        disp(['The file' f_name{control1} 'was skiped because no contains "raff" or "RAFF" in the name']);
                    end

                    if ~isempty(handles.SOURCE(control1))
                        % Enable fitting-button
                        set(handles.btn_fit_mt,'enable','on');
                        handles.SOURCE(control1).FTDATA=double(handles.SOURCE(control1).FTDATA);
                        % Disable output-options
                        set(handles.radio_saveopen_DICOM2,'enable','off');
                        set(handles.radio_saveopen_AEDES2,'enable','off');
                    end   
        end
        end
    % Case 2 the input files will be loaded RAFF (+Z) and (-Z) with 2 nifti files     
    case 2
            uiwait(msgbox('Please open multiples files with RAFF (+Z) and RAFF (-Z) in 2 nifti files from the next dialog!','Note','modal'));
            % Load RAFF data to handles structure FROM NIFTI SOURCE WITH TWO READS
            [f_name,f_path,f_index] = aedes_juigetfiles;
            n_files = size(f_name,2);
            control1= n_files;
            if (f_index ~= 0)
            % Go through all the file saving one by one in handles structure
            for control1 = 1 : control1
                % Skiping the files that has '_map' on name (output files).
                if(contains(f_name{control1},'_map')==1)
                disp(['The file ' f_name{control1} ' was skiped']);
                continue;
                % Saving the inversion files to list_minus (list of all
                % inversion files).
                elseif (contains(f_name{control1},'inv')== 1)
                    list_minus(i) = f_name(control1); 
                    i = i+1;
                    % Saving the plus files to list_minus (list of all
                    % inversion files).
                elseif(contains(f_name{control1},'plu')==1)
                    list_plus(p) = f_name(control1);
                    p = p+1;
                end
            end
                % Checking if the number of plus files and number of
                % minus files is diferent.
                if(i~=p)
                    uiwait(msgbox('Sorry, the number of files is wrong. Load RAFF (+Z) and RAFF (-Z) dataset','Report Window','modal'));
                    control3 = 2;
                    guidata(hObject,handles);
                end
                % Continous if (i=p)
                if(control3 ==1)
                    % Go through all the minus vector files
                    for i= 1:(size(list_minus,2))
                        % Delete the identifier '_inv' of the name without
                        % numeration(example:'_15_')
                        cell = strfind (list_minus{i},'_');
                        name_minus = erase(erase(list_minus{i},list_minus{i}(cell(1):cell(2))),'_inv');
                        % Go through all the plus vector files
                        for p = 1:(size(list_plus,2))
                            % Delete the identifier '_plu' of the name without
                            % numeration
                            cell_plus = strfind (list_plus{p},'_'); 
                            name_plus = erase(erase(list_plus{p},list_plus{p}(cell_plus(1):cell_plus(2))),'_plu');
                            % Check if the name of minus file is the same of the
                            % plus files and set datas on tmp1 and tmp2
                            if(contains(name_plus,name_minus)==1)
                                tmp1 = aedes_data_read([f_path{p},list_plus{p}]);
                                tmp2 = aedes_data_read([f_path{i},list_minus{i}]);
                                % Transform 2 files in one (tmp1)
                                if ~isempty(tmp2)
                                    tmp1.FTDATA=cat(ndims(tmp1.FTDATA),tmp1.FTDATA,tmp2.FTDATA);
                                end
                                % Load file to handles structure
                                handles.SOURCE(i)=tmp1;
                                if ~isempty(handles.SOURCE(i))
                                    % Enable fitting-button
                                    set(handles.btn_fit_mt,'enable','on');
                                    handles.SOURCE(i).FTDATA=double(handles.SOURCE(i).FTDATA);
                                    % Disable output-buttons
                                    set(handles.radio_saveopen_DICOM2,'enable','off');
                                    set(handles.radio_saveopen_AEDES2,'enable','off');
                                end       
                            end
                        end
                    end
                end
            end  
end

handles.LR=2;
handles.NiftiInput=1;
handles.SOURCETYPE='raff_map_';
% Stuff back to GUI for later use..
guidata(hObject,handles);






% --- Executes on button press in convert_2_1.
function convert_2_1_Callback(hObject, eventdata, handles)
% hObject    handle to convert_2_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This function is for users that want to join the (+Z) RAFF and (-Z) RAFF
% in one file, the output file is saving with the name of minus files with
% a sulfixe '_plus'
[f_name,f_path,f_index] = aedes_juigetfiles;
n_files = size(f_name,2);
control3 = 1;  
i = 1;
p = 1;
% Go through all the file selecting one by one to separate in 2 lists
for n_files = 1 : n_files
        if(contains(f_name{n_files},'_map')==1)
                disp(['The file ' f_name{n_files} ' was skiped']);
                continue;
        elseif (contains(f_name{n_files},'inv')== 1)
                list_minus(i) = f_name(n_files);
                i = i+1;
        elseif(contains(f_name{n_files},'plu')==1)
                list_plus(p) = f_name(n_files);
                p = p+1;
        end
 
end
        % Checking if the numbers of files in the lists is the same
        if(mod((i+p),2)==1)
            uiwait(msgbox('Sorry, the number of files is wrong. Load RAFF (+Z) and RAFF (-Z) dataset','Note!','modal'));
            control3 = 2;
            guidata(hObject,handles);
        end


        if(control3 ==1)
            for i= 1:(size(list_minus,2))
                % Delete the identifier '_inv' of the name without
                % numeration(example:'_15_')
                cell = strfind (list_minus{i},'_');
                name_minus = erase(erase(list_minus{i},list_minus{i}(cell(1):cell(2))),'_inv');
                % Go through all the plus vector files
                for p = 1:(size(list_plus,2))
                    % Delete the identifier '_plu' of the name without
                    % numeration
                    cell_plus = strfind (list_plus{p},'_'); 
                    name_plus = erase(erase(list_plus{p},list_plus{p}(cell_plus(1):cell_plus(2))),'_plu');
                    % Check if the name of minus file is the same of the
                    % plus files and set datas on tmp1 and tmp2
                    if(contains(name_minus,name_plus)==1)
                        tmp1 = aedes_data_read([f_path{p},list_plus{p}]); 
                        tmp2 = aedes_data_read([f_path{i},list_minus{i}]); 

                        if ~isempty(tmp2)
                            tmp1.FTDATA=cat(ndims(tmp1.FTDATA),tmp1.FTDATA,tmp2.FTDATA);
                        end
                        [fpath fname ext] = fileparts(list_minus{i}); %#ok<NCOMMA>
                        pre_fname{i} = [fname '_concatenated' ext]; %#ok<AGROW>
                        handles.SOURCE(i)=tmp1;
                    
                        if ~isempty(handles.SOURCE(i))
                            % Enable fitting-button
                            handles.SOURCE(i).FTDATA=double(handles.SOURCE(i).FTDATA);
                        end   

                    end

                end
            end
        end
        

        
% [fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
%                     'NIfTI Files - One File Format (*.nii)';...
%                     '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
%                     '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
%                     '*.*','All Files (*.*)'},'Save Data As',pre_fname);
for n_files= 1:size(list_plus,2)    
        % Then save MAP with given name!
        aedes_write_nifti(handles.SOURCE(n_files),fullfile(f_path{n_files},pre_fname{n_files}))
end
uiwait(msgbox('Files was sucessful converted','Note!','modal'));





%%%%%%%%%%%%%%%%%%%%%%%% FIT BUTTONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%
%
%


% --- Executes on button press in btn_fit_mt.
function btn_fit_mt_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fit_mt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if was loaded multiples files
if ((size(handles.SOURCE,2)) > 1);
    % Get the number of files
    n_files = size(handles.SOURCE,2);
    % Get the time sequence from user
    resp = aedes_inputdlg('Type tSL values','Input dialog',...
    {mat2str([])});    
    % Go through all the files
        for n_files = 1:n_files
            
            f_name = handles.SOURCE(n_files).HDR.fname;
            f_path = handles.SOURCE(n_files).HDR.fpath; 

            [dirPath, filename, ext] = fileparts(f_name);
            % Add the sulfixe '_map' to the result files
            result_name = [filename '_map' ext];
            result_path = f_path;
            Time_data = (size(handles.SOURCE(n_files).FTDATA))/2;
            Time_user = size(str2num(resp{1}));    
            % Check if the time sequence setted is the same of the file
            if (Time_data(4)) ~= (Time_user(2))
                disp(['The time sequence choose for ' f_name ' is different than the file'])
                continue;
            end
            
        
            if isfield(handles,'SOURCE') && ~isempty(handles.SOURCE(n_files))
            % Maybe there is enough data.. check fit type and noise masking and
            % output
            fit_mask=[get(handles.radiobutton11,'Value'),...
                get(handles.radiobutton12,'Value'),...
                get(handles.radiobutton13,'Value'),...
                get(handles.radiobutton14,'Value')];
            noise_mask=get(handles.checkbox2,'Value');
            output_mask=[get(handles.radio_saveopen_DICOM2,'Value'),...
                get(handles.radio_saveopen_AEDES2,'Value')    ];
    
            % check noise and make ROI if reasonable
            ROI=l_get_noise({handles.SOURCE(n_files)},noise_mask,handles.DEBUG);
            %keyboard;

            switch find(fit_mask)
                case 1  % 3-parameter non-linear            
                    handles.MAP=l_fit_MT({handles.SOURCE(n_files)},ROI,1,0,result_name,resp);

                case 2  % 4-parameter non-linear
                    handles.MAP=l_fit_MT({handles.SOURCE(n_files)},ROI,2,0,result_name,resp);

                case 3  % 3-parameter non-linear with T1
                    handles.MAP=l_fit_MT({handles.SOURCE(n_files)},ROI,3,0,result_name,resp);

                case 4  % 4-parameter non-linear with T1
                    if isfield(handles,'T1MAP')
                        % check & match sizes
                        [handles.T1MAP,stat]=l_checkT1size({handles.SOURCE(n_files)},{handles.T1MAP});
                        if stat==0
                            % something failed
                            disp('Couldn''t match T1 to MT/RAFF data. Check your data.');
                            return;
                        end
                        % T1 should now match MT/RAFF - fit map:

                        % ----------TODO-------------
                        % IMPLEMENT FITTING IN THE FUNCTION BELOW:



                        %handles.MAP=l_fit_MT({handles.SOURCE},ROI,4,{handles.T1MAP});




                        disp('T1-enabled fit not implemented yet!')
                        return;


                    else
                        disp('T1 map not loaded');
                    end


            end


    
            if ~isempty(handles.MAP)
                % Always save Nifti if enabled
                if get(handles.check_savenifti,'value')

                    % CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
                    % NOTE // THIS MAY CAUSE OTHER CONFUSION...
                    if ~isempty(handles.SOURCE(n_files)) && isfield(handles.SOURCE(n_files),'HDR') && isfield(handles.SOURCE(n_files).HDR,'fpath') && ~isempty(handles.SOURCE(n_files).HDR.fpath)
                        cd(handles.SOURCE(n_files).HDR.fpath);

                    end

                    % check if Nifti-input and replace header if yes:
                    if handles.NiftiInput
                        handles.MAP.DataFormat=handles.SOURCE.DataFormat;
                        handles.MAP.HDR=handles.SOURCE.HDR;
                    end

                    % figure out a pre-suggestion for filename
                    nn=1;
                    pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);

                    while exist(pre_fname)
                        nn=nn+1;
                        pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);
                    end

                    % pre-suggestion ready, give dialog with it:
    %                 [fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
    %                     'NIfTI Files - One File Format (*.nii)';...
    %                     '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
    %                     '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
    %                     '*.*','All Files (*.*)'},'Save Data As',pre_fname);

                    if ~isequal(f_name,0) || isequal(f_path,0)

                        % Then save MAP with given name!
                        aedes_write_nifti(handles.MAP,fullfile(result_path,result_name));
                    end

                end % if save-nifti


    %             switch find(output_mask)
    %                 case 1 % export as DICOM
    %                     l_export_dicom({handles.MAP});
    %                 case 2 % open in AEDES
    %                     aedes(handles.MAP);
    %             end

                % Stuff back to GUI for later use..
                guidata(hObject,handles);

            end
            % keyboard;

        else
            disp('Please load data first!');
            return;
        end

        end 

else
    if isfield(handles,'SOURCE') && ~isempty(handles.SOURCE)
        % Maybe there is enough data.. check fit type and noise masking and
        % output
        fit_mask=[get(handles.radiobutton11,'Value'),...
        get(handles.radiobutton12,'Value'),...
        get(handles.radiobutton13,'Value'),...
        get(handles.radiobutton14,'Value')];
        noise_mask=get(handles.checkbox2,'Value');
        output_mask=[get(handles.radio_saveopen_DICOM2,'Value'),...
        get(handles.radio_saveopen_AEDES2,'Value')    ];

        % check noise and make ROI if reasonable
        ROI=l_get_noise({handles.SOURCE},noise_mask,handles.DEBUG);
        %keyboard;

        switch find(fit_mask)
            case 1  % 3-parameter non-linear            
                handles.MAP=l_fit_MT({handles.SOURCE},ROI,1);

            case 2  % 4-parameter non-linear
                handles.MAP=l_fit_MT({handles.SOURCE},ROI,2);

            case 3  % 3-parameter non-linear with T1
                handles.MAP=l_fit_MT({handles.SOURCE},ROI,3);

            case 4  % 4-parameter non-linear with T1
                if isfield(handles,'T1MAP')
                        % check & match sizes
                    [handles.T1MAP,stat]=l_checkT1size({handles.SOURCE},{handles.T1MAP});
                    if stat==0
                        % something failed
                        disp('Couldn''t match T1 to MT/RAFF data. Check your data.');
                        return;
                    end
                    % T1 should now match MT/RAFF - fit map:

                    % ----------TODO-------------
                    % IMPLEMENT FITTING IN THE FUNCTION BELOW:
                    %handles.MAP=l_fit_MT({handles.SOURCE},ROI,4,{handles.T1MAP});

                    disp('T1-enabled fit not implemented yet!')
                    return;


                else
                    disp('T1 map not loaded');
                end


        end
        if ~isempty(handles.MAP)
            % Always save Nifti if enabled
            if get(handles.check_savenifti,'value')

            % CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
            % NOTE // THIS MAY CAUSE OTHER CONFUSION...
                if ~isempty(handles.SOURCE) && isfield(handles.SOURCE,'HDR') && isfield(handles.SOURCE.HDR,'fpath') && ~isempty(handles.SOURCE.HDR.fpath)
                    cd(handles.SOURCE.HDR.fpath);
                end

                % check if Nifti-input and replace header if yes:
                if handles.NiftiInput
                    handles.MAP.DataFormat=handles.SOURCE.DataFormat;
                    handles.MAP.HDR=handles.SOURCE.HDR;
                end

                % figure out a pre-suggestion for filename
                nn=1;
                pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);

                while exist(pre_fname)
                    nn=nn+1;
                    pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);
                end

                % pre-suggestion ready, give dialog with it:
                [fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
                'NIfTI Files - One File Format (*.nii)';...
                '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
                '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
                '*.*','All Files (*.*)'},'Save Data As',pre_fname);

                if ~isequal(fname,0) || isequal(fpath,0)
                    % Then save MAP with given name!
                    aedes_write_nifti(handles.MAP,fullfile(fpath,fname));
                end

            end % if save-nifti


            switch find(output_mask)
                case 1 % export as DICOM
                        l_export_dicom({handles.MAP});
                case 2 % open in AEDES
                        aedes(handles.MAP);
            end

                % Stuff back to GUI for later use..
                guidata(hObject,handles);

        end
            % keyboard;

        else
            disp('Please load data first!');
            return;
    end
end







% --- Executes on button press in btn_fit_t1r.
function btn_fit_t1r_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fit_t1r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Check if was loaded multiples files
if size(handles.SOURCE,2) > 1
    % Get the number of files       
    n_files = size(handles.SOURCE,2);
    % Get the time sequence from user
    resp = aedes_inputdlg('Type tSL values','Input dialog',...
    {mat2str([])});

    % Go through all the files
    for n_files = 1 : n_files  
    
        f_name = handles.SOURCE(n_files).HDR.fname;
        f_path = handles.SOURCE(n_files).HDR.fpath; 

        T1r_name = contains(f_name,'T1r');
        T2r_name = contains(f_name,'T2r');
        map_name = contains(f_name,'_map');
        % Verify if the file has 'T1r' or 'T2r' in name, if has '_map'
        % skipped
        if (T1r_name == 0) && (T2r_name == 0) || (map_name == 1)
            disp(['The ' f_name ' was skipped']);
            continue;
        end
        
        % Check if the time sequence setted is the same of the file
        Time_data = size(handles.SOURCE(n_files).FTDATA);
        Time_user = size(str2num(resp{1}));    
    
        if (Time_data(4)) ~= (Time_user(2))
            disp(['The sequences choose for ' f_name ' is different than the file'])
            continue;
        end


        %if ~isempty(handles.SOURCE)
        %   % Enable fitting-button
        %   set(handles.btn_fit_t1r,'enable','on');
        %   handles.SOURCE.FTDATA=double(handles.SOURCE.FTDATA);
        %end

        % set a tag for identifying MT/RAFF vs. T12rho
        handles.LR=1;
        handles.NiftiInput=1;
        handles.SOURCETYPE='t1rho_map_';
        
        % Add the sulfixe '_map' to the result files
        [dirPath, filename, ext] = fileparts(f_name);
        result_name = [filename '_map' ext];
        result_path = f_path;
 

        if isfield(handles,'SOURCE') && ~isempty(handles.SOURCE(n_files))
            %Maybe there is enough data.. check fit type and noise masking and output
            fit_mask=[get(handles.radiobutton1,'Value'),...
                get(handles.radiobutton2,'Value'),...
                get(handles.radiobutton3,'Value'),...
                get(handles.radiobutton4,'Value')];
            noise_mask=get(handles.checkbox1,'Value');
            output_mask=[get(handles.radio_saveopen_DICOM1,'Value'),...
            get(handles.radio_saveopen_AEDES1,'Value')    ];
    
            %  check noise and make ROI if reasonable
            ROI=l_get_noise({handles.SOURCE(n_files)},noise_mask,handles.DEBUG);
            % keyboard;
   


            switch find(fit_mask)
                case 1  % 2-parameter linear            
                    handles.MAP=l_fit_T1r({handles.SOURCE(n_files)},ROI,1,0,result_name,resp);
            
                case 2  % 2-parameter non-linear
                    handles.MAP=l_fit_T1r({handles.SOURCE(n_files)},ROI,2,0,result_name,resp);

                case 3  % 3-parameter non-linear
                    handles.MAP=l_fit_T1r({handles.SOURCE(n_files)},ROI,3,0,result_name,resp);

                case 4
                    if isfield(handles,'T1MAP')
                        % check & match sizes
                        [handles.T1MAP,stat]=l_checkT1size({handles.SOURCE(n_files)},{handles.T1MAP});
                        if stat==0
                            % something failed
                            disp('Couldn''t match T1 to T1/2rho data. Check your data.');
                            return;
                        end
                    % T1 should now match T1r - fit map:
                
                    % ----------TODO-------------
                    % IMPLEMENT FITTING IN THE FUNCTION BELOW:

                        handles.MAP=l_fit_T1r({handles.SOURCE(n_files)},ROI,4,{handles.T1MAP});                
                    else
                    M     disp('T1 map not loaded');
                    end

            end
    
            if ~isempty(handles.MAP)   % possible to be empty if manual values requested but none given...
                % Always save Nifti if enabled
                if get(handles.check_savenifti,'value')
            
                    % CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
                    % NOTE // THIS MAY CAUSE OTHER CONFUSION...
                    if ~isempty(handles.SOURCE(n_files)) && isfield(handles.SOURCE(n_files),'HDR') && isfield(handles.SOURCE(n_files).HDR,'fpath') && ~isempty(handles.SOURCE(n_files).HDR.fpath)
                        cd(handles.SOURCE(n_files).HDR.fpath);
                    end
            
                    % check if Nifti-input and replace header if yes:
                    if handles.NiftiInput
                        handles.MAP.DataFormat=handles.SOURCE.DataFormat;
                        handles.MAP.HDR=handles.SOURCE.HDR;
                    end
            
                    % figure out a pre-suggestion for filename
                    nn=1;
                    pre_fname=sprintf('%s_%02d.nii',result_name,nn);
            
                    while exist(pre_fname)
                        nn=nn+1;
                        pre_fname=sprintf('%s%_02d.nii',result_name,nn);
                    end
            
                    % pre-suggestion ready, give dialog with it:
                    %[fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
                    %    'NIfTI Files - One File Format (*.nii)';...
                    %   '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
                    %    '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
                    %    '*.*','All Files (*.*)'},'Save Data As',pre_fname);
            
                    if ~isequal(f_name,0) || isequal(f_path,0)
                        % Then save MAP with given name!
                        aedes_write_nifti(handles.MAP,fullfile(result_path,result_name));
                    end
            
                end % if save-nifti

        
                switch find(output_mask)
                    case 1 % export as DICOM
                        %l_export_dicom({handles.MAP});
                    case 2 % open in AEDES
                        %aedes(handles.MAP);
                end
     
                % Stuff back to GUI for later use..
            
        
            end
    
        % keyboard;
    
        else
            disp('Please load data first!');
            return;
        end




    end
    % Stuff back to GUI for later use..
    guidata(hObject,handles);
    
else
    if isfield(handles,'SOURCE') && ~isempty(handles.SOURCE)
        % Maybe there is enough data.. check fit type and noise masking and
        % output
        fit_mask=[get(handles.radiobutton1,'Value'),...
            get(handles.radiobutton2,'Value'),...
            get(handles.radiobutton3,'Value'),...
            get(handles.radiobutton4,'Value')];
        noise_mask=get(handles.checkbox1,'Value');
        output_mask=[get(handles.radio_saveopen_DICOM1,'Value'),...
            get(handles.radio_saveopen_AEDES1,'Value')    ];
        
        % check noise and make ROI if reasonable
        ROI=l_get_noise({handles.SOURCE},noise_mask,handles.DEBUG);
        %keyboard;
    
        switch find(fit_mask)
            case 1  % 2-parameter linear            
                handles.MAP=l_fit_T1r({handles.SOURCE},ROI,1);
                
            case 2  % 2-parameter non-linear
                handles.MAP=l_fit_T1r({handles.SOURCE},ROI,2);

            case 3  % 3-parameter non-linear
                handles.MAP=l_fit_T1r({handles.SOURCE},ROI,3);

            case 4
                if isfield(handles,'T1MAP')
                    % check & match sizes
                    [handles.T1MAP,stat]=l_checkT1size({handles.SOURCE},{handles.T1MAP});
                    if stat==0
                        % something failed
                        disp('Couldn''t match T1 to T1/2rho data. Check your data.');
                        return;
                    %[fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
                    %'NIfTI Files - One File Format (*.nii)';...
                    %'*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
                    %'*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
                    %'*.*','All Files (*.*)'},'Save Data As',pre_fname);
                    end
                    % T1 should now match T1r - fit map:
                
                    % ----------TODO-------------
                    % IMPLEMENT FITTING IN THE FUNCTION BELOW:

                    handles.MAP=l_fit_T1r({handles.SOURCE},ROI,4,{handles.T1MAP});           
                
                else
                    disp('T1 map not loaded');
                end     
        end
    
        if ~isempty(handles.MAP)   % possible to be empty if manual values requested but none given...
            % Always save Nifti if enabled
            if get(handles.check_savenifti,'value')
            
                % CHANGE DIRECTORY TO THE FIRST FILE DIRECTORY..
                % NOTE // THIS MAY CAUSE OTHER CONFUSION...
                if ~isempty(handles.SOURCE) && isfield(handles.SOURCE,'HDR') && isfield(handles.SOURCE.HDR,'fpath') && ~isempty(handles.SOURCE.HDR.fpath)
                    cd(handles.SOURCE.HDR.fpath);
                
                end
            
                % check if Nifti-input and replace header if yes:
                if handles.NiftiInput
                    handles.MAP.DataFormat=handles.SOURCE.DataFormat;
                    handles.MAP.HDR=handles.SOURCE.HDR;
                end
            
                % figure out a pre-suggestion for filename
                nn=1;
                pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);
            
                while exist(pre_fname)
                    nn=nn+1;
                    pre_fname=sprintf('%s%02d.nii',handles.SOURCETYPE,nn);
                end
            
                % pre-suggestion ready, give dialog with it:
                [fname,fpath,findex]=uiputfile({'*.nii;*.NII',...
                'NIfTI Files - One File Format (*.nii)';...
                '*.hdr;*.HDR','NIfTI Files - Two File Format (*.hdr)';...
                '*.hdr;*.HDR','Analyze 7.5 Files (*.hdr)';...
                '*.*','All Files (*.*)'},'Save Data As',pre_fname);           
                if ~isequal(fname,0) || isequal(fpath,0)
                
                    % Then save MAP with given name!
                    aedes_write_nifti(handles.MAP,fullfile(fpath,fname));
                end
            
            end % if save-nifti

        
            switch find(output_mask)
                case 1 % export as DICOM
                    l_export_dicom({handles.MAP});
                case 2 % open in AEDES
                    aedes(handles.MAP);
            end
        
        
            % Stuff back to GUI for later use..
            guidata(hObject,handles);
            
        
        end
    
    % keyboard;
    
    else
        disp('Please load data first!');
        return;
    end
    
end










% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2

%keyboard;
handles.DEBUG=get(handles.togglebutton2,'value');

if handles.DEBUG
    set(handles.togglebutton2,'foregroundColor',[1 1 1]);
    set(handles.togglebutton2,'backgroundColor',[1 0 0]);
    set(handles.togglebutton2,'string','debug ON');
else
    set(handles.togglebutton2,'foregroundColor',[0 0 0]);
    set(handles.togglebutton2,'backgroundColor',[0.941176 0.941176 0.941176]);
    set(handles.togglebutton2,'string','debug off');
end

% As it is run, stuff handles - structure to workspace..
assignin('base','handles',handles);

guidata(hObject,handles);





% Check ROI MASKS
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'SOURCE') && ~isempty(handles.SOURCE)
    % Maybe there is enough data.. check fit type and noise masking and
    % output
    
    if handles.LR==1
        fit_mask=[get(handles.radiobutton1,'Value'),...
            get(handles.radiobutton2,'Value'),...
            get(handles.radiobutton3,'Value'),...
            get(handles.radiobutton4,'Value')];
        noise_mask=get(handles.checkbox1,'Value');
        output_mask=[get(handles.radio_saveopen_DICOM1,'Value'),...
            get(handles.radio_saveopen_AEDES1,'Value')    ];
    
        % check noise and make ROI if reasonable
        ROI=l_get_noise({handles.SOURCE},noise_mask,handles.DEBUG);
    elseif handles.LR==2
        fit_mask=[get(handles.radiobutton11,'Value'),...
            get(handles.radiobutton12,'Value'),...
            get(handles.radiobutton13,'Value'),...
            get(handles.radiobutton14,'Value')];
        noise_mask=get(handles.checkbox2,'Value');
        output_mask=[get(handles.radio_saveopen_DICOM2,'Value'),...
            get(handles.radio_saveopen_AEDES2,'Value')    ];
        
        %keyboard;
        % check noise and make ROI if reasonable
        ROI=l_get_noise({handles.SOURCE},noise_mask,handles.DEBUG);
    else
        % something wrong..
        warndlg('Couldn''t do something. Try restarting program.','Problem');
    end
    
      %  keyboard;
    
    
    % Save ROIs first:
    % Prompt for file name
    %[fname,fpath,findex]=uiputfile({'*.roi' ,...
    %    'ROI-Files (*.roi)';...
    %    '*.*','All Files (*.*)'},...
    %    'Save ROI-file',['guessed_masks.roi']);
    %if isequal(fname,0) || isequal(fpath,0)
    %    return
    %end
    %[fp,fn,fe]=fileparts(fname);
    %if isempty(fe)
    %    fe='.roi';
    %end
    %save(fullfile(fpath,[fn fe]),'ROI','-mat');
    
    DATA=handles.SOURCE;
    DATA.FTDATA=DATA.FTDATA(:,:,:,1);
       
    aedes(DATA,ROI);
    
else
    % no source data.. nothing to do
     warndlg('No source data, can'' work.','No data');
end










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   LOCAL FUNCTIONS
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%






function [f,stat]=l_checkT1size(T1r,T1map)
if ~(ndims(T1r{1}.FTDATA(:,:,:,1))==ndims(T1map{1}.FTDATA))
    disp('Dimensions not matching, aborting.');
    f=[];stat=0;
    return;
end

if size(T1r{1}.FTDATA,3)>size(T1map{1}.FTDATA,3)
    disp('T1map too small to use!');
    f=[];stat=0;
    return;
end
    
% what else could fail? A number of things clearly... but let's try now.
if sum(abs(size(T1r{1}.FTDATA(:,:,:,1))-size(T1map{1}.FTDATA)))==0
    % easiest case: size equals size done already!
    f=T1map{1};
    stat=1;
    return
else
    % some sizes are different. Let's resize
    
    
    % -------FIXME--------- 
    % STARTING ASSUMPTION: T1 has more slices
    % T1 may have larger matrix size
    
    % slice indices for T1:
    ss=[];
    for ii=1:length(T1map{1}.HDR.dicominfo)
        ss(ii)=T1map{1}.HDR.dicominfo(ii).SliceLocation;
    end
    % slice indices for T1rho
    ssT1r=[];
    for ii=1:size(T1r{1}.FTDATA,3)
        ssT1r(ii)=T1r{1}.HDR.dicominfo(ii).SliceLocation;
    end
    % starting index for minimum difference between the two:
    sti=[];
    for ii=1:length(ss)-length(ssT1r)
        sti(ii)=sum(ss(ii:ii+length(ssT1r)-1)-ssT1r);
    end
    
    % And the minimum should give the starting index:
    [dumm,sti_T1]=min(abs(sti));
    
    % Remove extra slices from T1 map:
    T1map{1}.FTDATA=T1map{1}.FTDATA(:,:,sti_T1:sti_T1-1+length(ssT1r));
    
    % Re-check the dimensions..
    if sum(abs(size(T1r{1}.FTDATA(:,:,:,1))-size(T1map{1}.FTDATA)))==0
        % easiest case: size equals size done already!
        f=T1map{1};
        stat=1;
        return;
    else
        % Ok - NOW resize
        tmp=[];
        for ii=1:length(ssT1r)
            tmp(:,:,ii)=imresize(T1map{1}.FTDATA(:,:,ii),[size(T1r{1}.FTDATA,1) size(T1r{1}.FTDATA,2)]);
        end
        T1map{1}.FTDATA=tmp;
        f=T1map{1};
        stat=1;
        return;
    end
    

end





function ROI=l_get_noise(DATA,noise_mask,ddebug)

%keyboard;

if nargin<3
    ddebug=0;
end

ROI=[];
if noise_mask==0
        return;
end

   
% Figure out 3-D (2-D multislice) vs. 2-D
% number of slices - from d_loader, the volumes contain TE or TR /whatnot
% and the dim(3) contains slices, IF there are 4 dimensions to the data
if (ndims(DATA{1}.FTDATA))==4
    if size(DATA{1}.FTDATA,4)==2
        % this is most likely single-slice data
        ns=1;
        
        % RESHAPE DATA TO THE SAME FORMAT AS MULTISLICE:
        % [x, y, slices, tSL]
        szm=size(DATA{1}.FTDATA);
        DATA{1}.FTDATA=reshape(DATA{1}.FTDATA,[szm(1),szm(2),1,szm(3)*szm(4)]);
    end
    % Now dim3 should contain correct number of slices
    ns=size(DATA{1}.FTDATA,3);
else
    ns=1;
end

if ns>1
    ns=2;
end



    ROI(1).label='ROI 1';
    ROI(2).label='bg';
    ROI(1).color=[255 0 0];
    ROI(2).color=[0 255 0];
    % create template ROI using first volume:
    ROI(1).voxels{1}=logical(ones(size(DATA{1}.FTDATA(:,:,:,1))));
    ROI(2).voxels{1}=logical(ones(size(DATA{1}.FTDATA(:,:,:,1))));
    
    % this will generate partially aedes-like ROI
    switch ns
        case 1  % 2-D + decay, to be implemented
            % go through slices and generate ROIs for each
            for ii=1:size(ROI(1).voxels{1},3)
                tmp=double(DATA{1}.FTDATA(:,:,ii));
                [aa,bb]=hist(tmp(:),60); % take 60-bin histogram
                [zc,izc]=l_zerocross(diff(aa)); % find zerocrossings of the 1st derivate
                % assume noise has the first peak and data the second ->
                % something b/w first point and first dip likely a good
                % threshold value
                %keyboard;
                th=bb(floor((1+izc(1))/2));
                mask=ones(size(tmp));
                mask(tmp<th)=0; % this is the image mask
                ROI(2).voxels{1}(:,:,ii)=logical(1-mask); % noise/bg mask
                
                if ddebug
                    mask=zeros(size(mask));
                    sz1=floor(size(mask,1)/2);sz2=floor(size(mask,2)/2);
                    ms=5;
                    mask(sz1-ms:sz1+ms-1,sz2-ms:sz2+ms-1)=ones(ms*2,ms*2);
                end
                
                ROI(1).voxels{1}(:,:,ii)=logical(mask); % fit mask
            end
            
        case 2
            % 3-D + dacey.  go through the first volume and generate ROI
            
            % TAKE THE FIRST VOLUME!
            tmp=double(DATA{1}.FTDATA(:,:,:,1));
            [aa,bb]=hist(tmp(tmp>=0),60); % take 60-bin histogram % DM modify ---was  [aa,bb]=hist(tmp(:),60);
            [zc,izc]=zerocross(diff(aa)); % find zerocrossings of the 1st derivate
            % assume noise has the first peak and data the second ->
            % something b/w first point and first dip likely a good
            % threshold value
            % keyboard;
            %th=bb(floor((1+izc(1))/2));
            th=1;
            mask=ones(size(tmp));
            mask(tmp<th)=0; % this is the image mask
            ROI(2).voxels{1}(:,:,:)=logical(1-mask); % noise/bg mask
            
            if ddebug
                mask=zeros(size(mask));
                sz1=floor(size(mask,1)/2);sz2=floor(size(mask,2)/2);
                ms=5;
                mask(sz1-ms:sz1+ms-1,sz2-ms:sz2+ms-1)=ones(ms*2,ms*2);
            end
            
            ROI(1).voxels{1}(:,:,:)=logical(mask); % fit mask
         
            
            
        otherwise
            % go through slices and generate ROIs for each
            for ii=1:size(ROI(1).voxels{1},3)
                tmp=double(DATA{1}.FTDATA(:,:,ii,1));                
                [aa,bb]=hist(tmp(:),60); % take 60-bin histogram
                [zc,izc]=l_zerocross(diff(aa)); % find zerocrossings of the 1st derivate
                % assume noise has the first peak and data the second ->
                % something b/w first point and first dip likely a good
                % threshold value
                %keyboard;
                %th=bb(floor((1+izc(1))/2));
                th=1;
                mask=ones(size(tmp));
                mask(tmp<th)=0; % this is the image mask
                ROI(2).voxels{1}(:,:,ii)=logical(1-mask); % noise/bg mask
                
                if dDebug
                    mask=zeros(size(mask));
                    sz1=floor(size(mask,1)/2);sz2=floor(size(mask,2)/2);
                    ms=5;
                    mask(sz1-ms:sz1+ms-1,sz2-ms:sz2+ms-1)=ones(ms*2,ms*2);
                end
                
                ROI(1).voxels{1}(:,:,ii)=logical(mask); % fit mask
            end
    end

    
    
    
function [f,ndx] = l_zerocross(vector)
%ZEROCROSS Finds number of zerocrosses
%    Zerocross simply reports the number of times
%    the input vector crosses the zero boundary. If two
%    outputs defined, also returns indices
%    (the latter index of crossing).
%
%    Usage
%      [NUMBER,INDEX]=ZEROCROSS(X);

% (c) mikko nissi 2012, <nissi@cmrr.umn.edu>

ndx=find(abs(diff(sign(vector))));
if~isempty(ndx)
    ndx=ndx+1;
    f=length(ndx);
else
    return
end



function l_export_dicom(DATA)
% This Aedes plugin exports FTDATA as dicom file

% This function is written for Aedes
%
% Updates
% 6/3/2013, Mikko Nissi
%    - Add search for data range (if T1rho map, etc) and re-scale!
%    - Add dicom-header-stuffing if originating from Dicom..
%
% copyright (c) 2011 Mikko Nissi <mikko.nissi@iki.fi>
%
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%keyboard;

% Try to find procpar if not included already..!
try
    % dicom-flag
    wasdicom=0;
    
    
    if isfield(DATA{1},'PROCPAR')
        % Ok, no worries..
        disp('Appears to be Varian data');
        p=DATA{1}.PROCPAR;
    elseif isfield(DATA{1},'HDR')&& (exist(fullfile(DATA{1}.HDR.fpath,'procpar'))==2)
        % try to load procpar
        p=aedes_readprocpar(fullfile(DATA{1}.HDR.fpath,'procpar'));
    elseif isfield(DATA{1},'HDR') && isfield(DATA{1}.HDR,'dicominfo')
        disp('Appears to be Dicom data');
        % File has dicom origin, use all of that!
        % set flag
        wasdicom=1;
        isparam=0;
        if isfield(DATA{1},'Param') && isfield(DATA{1}.Param,'Header')
            % seems to be a parametric map - set another flag
            disp('Appears to be parametric map also');
            isparam=1;
        end
        
        
    else
        p.name={'Unknown'};
        p.seqfil={'Unknown sequence'};
        p.comment={'no comment'};
        p.np=2;p.lro=.1; % should yield 1
        p.nv=1;p.lpe=.1; % should yield 1
        p.dimX={'lpe'};
        p.dimY={'lro'};
        p.operator_={'John Doe'};
        p.date={'2011'};
        p.reffrq='400.1';
        p.tn={'H1'};
        p.rfcoil={'rfcoil'};
        p.console={'nmr'};
        p.thk=1; % default to 1 mm
        p.acqdim=2;
    end
    
    % Collect data from procpar
    if wasdicom==0
        % Non-dicom origin, use guessed / varian-derived stuff for dicom
        patientsname=p.name{1};
        seriesdescription=sprintf('%s %s',p.seqfil{1},p.comment{1});
        performingphysiciansname=p.operator_{1};
        imagingfrequency=p.reffrq;
        imagednucleus=p.tn{1};
        studydate=p.date{1};
        softwareversions=p.console{1};
        transmitcoilname=p.rfcoil{1};
        slicethickness=p.thk;
        %rows=p.nv;
        %cols=p.np/2;
        orient=[0;1;0;0;1;0]; % some random fixed orientation parameter...
        patpos=[0;0;0];
        
        if strcmpi(p.dimX,'lpe')
            pixelspacing=[p.lpe/p.nv*10 p.lro/(p.np/2)*10];
        elseif strcmpi(p.dimX,'lro')
            pixelspacing=[p.lro/(p.np/2)*10 p.lpe/p.nv*10];
        else
            pixelspacing=[1 1];
        end
        if p.acqdim==3
            slicethickness=p.lpe2/(p.nv2)*10;
            mracqtype='3D';
        else
            mracqtype='2D';
        end
        
        
        Dat=double(DATA{1}.FTDATA);

        % One final check for parametric data
        % If there''s "param" -field, do not scale, but cut the values to dicom range
        if isfield(DATA{1},'Param')
            % seems to be a parametric map - set another flag
            Dat(Dat>4095)=4095;
            Dat(Dat<0)=0;
            Dat=uint16(Dat);
        else
                   
            % Scale data (linear!!)
            Dat=Dat-min(Dat(:));
            Dat=Dat./max(Dat(:));
            Dat=uint16(Dat*4095);
        end
       
    else
        % Set parameters from dicom        
        % See if there's a scale factor..
        rescale=0;
        
        if isfield(DATA{1},'Param') && isfield(DATA{1}.Param,'Range')
            Dat=double(DATA{1}.FTDATA);
            lims=DATA{1}.Param(1).Range;
            Dat(Dat<lims(1))=lims(1);
            Dat(Dat>lims(2))=lims(2);
            
         
            % Scale data (linear!!)
            Dat=Dat-lims(1);
            Dat=Dat./lims(2);
            Dat=uint16(Dat*4095);
            
            Rescale_intercept=lims(1);
            Rescale_slope=lims(2)/4095;
            
            rescale=1;
        else
            % if it wasn't parametric, keep whatever scale there was in dicom!
            Dat=DATA{1}.FTDATA;
            
            if ~isinteger(Dat)
                % Convert to uint16 if needed
                disp('Dicom data, but values not integer');
                if max(Dat(:))>4095
                    disp(' \Outside Dicom range, re-scaling to 14 bit and converting')
                    % Rescale if out of bounds!
                    % Scale data (linear!!)
                    Dat=Dat-min(Dat(:));
                    Dat=Dat./max(Dat(:));
                    Dat=uint16(Dat*4095);
                else
                    disp(' \Converting to uint16');
                    Dat=uint16(Dat);
                end
            end
            
        end
       
    end
    
    
    % Prompt for file name
    [fname,fpath,findex]=uiputfile({'*.dcm' ,...
        'DICOM-Files (*.dcm)';...
        '*.*','All Files (*.*)'},...
        'Save DICOM-file',[DATA{1}.HDR.fpath, ...
        'dicomdata']);
    if isequal(fname,0) || isequal(fpath,0)
        return
    end
    
    tic
    
    [fp,fn,fe]=fileparts(fname);
    if isempty(fe)
        fe='.dcm';
    end
    
    %sprintf('%s_%03d%s',fullfile(fpath,fn),1,fe)
    
    % Save file 1 and then grab as template.
    if wasdicom==0
    
    dicomwrite(Dat(:,:,1),sprintf('%s_%03d%s',fullfile(fpath,fn),1,fe),...
        'Modality','MR',...
        'PatientName',patientsname,...
        'SeriesDescription',seriesdescription,...
        'PerformingPhysiciansName',performingphysiciansname,...
        'ImagingFrequency',imagingfrequency,...
        'ImagedNucleus',imagednucleus,...
        'SoftwareVersions',softwareversions,...
        'TransmitCoilName',transmitcoilname,...
        'PixelSpacing',pixelspacing,...
        'SpacingBetweenSlices',slicethickness,...
        'SeriesNumber',1,...
        'SliceThickness',slicethickness,...
        'MRAcquisitionType',mracqtype,...
        'PatientPosition','FFS',...
        'ImageOrientationPatient',orient,...
        'ImagePositionPatient',patpos);
    
    template=dicominfo(sprintf('%s_%03d%s',fullfile(fpath,fn),1,fe));
    
    end
    
    nslices=size(Dat,3);
    
    
    %keyboard;
    
    h=waitbar(0,'Exporting image');
    drawnow;
    
    for ii=1:nslices
        if wasdicom
            % use previous dicom-data
            if isparam
                % parametric map, use header from
                template=DATA{1}.Param(ii).Header(1);
                template.SeriesDescription=DATA{1}.Param(ii).Type;
                template.ProtocolName=DATA{1}.Param(ii).Type;
                % clear up some seriesnumber to separate from other data...
                template.SeriesNumber=template.SeriesNumber+100;  % ADD RANDOM LARGE NUMBER!!
                if rescale
                    % Add scaling factors!
                    template.RescaleIntercept=Rescale_intercept;
                    template.RescaleSlope=Rescale_slope;
                end
                
            else
                % just a dicom, use header from:
                template=DATA{1}.HDR.dicominfo(ii);
                template.SeriesNumber=template.SeriesNumber+100;  % ADD "RANDOM" LARGE NUMBER!!
            end
            %keyboard;
            dicomwrite(Dat(:,:,ii),sprintf('%s_%03d%s',fullfile(fpath,fn),ii,fe),template,'CreateMode','copy');
            waitbar(ii/nslices,h,sprintf('Exporting image %03d/%03d',ii,nslices));
            
        else
            % use the pre-save-template 
            template.SliceLocation=ii-nslices/2; % set zero somewhere in the middle..
            template.Filename=sprintf('%s_%03d%s',fullfile(fpath,fn),ii,fe);
            template.InstanceNumber=ii;
            dicomwrite(Dat(:,:,ii),sprintf('%s_%03d%s',fullfile(fpath,fn),ii,fe),template);
            waitbar(ii/nslices,h,sprintf('Exporting image %03d/%03d',ii,nslices));
        end
    end
    close(h);
    
    
    
catch
    errordlg({'Could not export files. The following error was returned',...
        '',lasterr},'modal')
end


toc




    
    
function f=l_fit_T1r(DATA,ROI,fit_type,T1DATA,savefname,resp)
% This Aedes plugin calculates t1/2rho map from Siemens-t1/2r-prep data

% This function is written for Aedes
%
% copyright (c) 2012 Mikko Nissi <nissi@cmrr.umn.edu>
%
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

% check input for T1DATA
if (nargin == 6)
    clearvars T1DATA;
    T1DATA{1}.FTDATA=zeros(size(DATA{1}.FTDATA(:,:,:,1)));
end
if nargin<4
    T1DATA{1}.FTDATA=zeros(size(DATA{1}.FTDATA(:,:,:,1)));   
end



if isfield(DATA{1}.HDR,'dicominfo')
    % DICOM FOUND, USE IT!


    
    % Grab dicom stuff
    p=DATA{1}.HDR.dicominfo;
    
    % Dig TR from PROCPAR (in seconds)
    TR = p(1).RepetitionTime;
    
    % dig sWiPMemBlock:
    [fl,fd,strrs]=l_new_parse_swip(DATA);
    
    % WIP-data, naming follows sequence internal variable names
    WIP_RectPulseDur = fl(1);
    WIP_NumPrepPulses= fl(2);
    WIP_PrepPulseDur= fl(3);
    WIP_T2rhoPulseDur= fl(4);
    WIP_InversionCheckBox= fl(5);
    WIP_InvPulseDur= fl(6);
    WIP_SpoilerMom= fl(7);
    WIP_MTOffsetFreq= fl(8);
    WIP_T2rhoCheckBox= fl(9);
    WIP_ArrayCheckBox= fl(10);
    WIP_InternalTR= fl(11);
    WIP_PrepPulseShape= fl(12);
    WIP_RAFFPulseDur= fl(13);
    WIP_NumRectPulses= fl(14);
    %WIP_MinIPD= fl(15);
    
    dWIP_PrepVoltageAdj = fd(1);
    dWIP_RAFFVoltageAdj= fd(2);
    dWIP_RectVoltageAdj= fd(3);
    dWIP_T2rVoltageAdj= fd(4);
    dWIP_InvVoltageAd= fd(5);
    
    
    pw=WIP_PrepPulseDur;  % us, T2r PREP pulse duration
    pw=pw/1000; % ms
    
    % For now, the array of pulses is as follows:
    %%% PLEASE NOTE ---- THIS DEPENDS ON THE NBPP parameter in the seq code
    if (WIP_ArrayCheckBox==1)
        n_mlev=[0:4:WIP_NumPrepPulses];
    end
    
    % create tSL-vector
    tSL=pw.*n_mlev;

    
else
    % NO DICOM.. display dialog
    if nargin ==3
    resp = aedes_inputdlg('Type tSL values','Input dialog',...
        {mat2str([])});
    end
    if isempty(resp)
        % NO RESPONSE OR SOME OTHER ISSUES.. EXIT
        f=[];
        return
    else
        resp=resp{1};
        fit_vals = str2num(resp);
        tSL=fit_vals;
    end
end

    
% averages: this has no meaning here..
nex=1;

% number of slices - from d_loader, the volumes contain TE or TR /whatnot
% and the dim(3) contains slices, IF there are 4 dimensions to the data
if (ndims(DATA{1}.FTDATA))==4
    ns=size(DATA{1}.FTDATA,3);
else
    ns=1;
end



%keyboard; %% find background ROI
% if masknoise, then signals below noise-threshold are not considered for
% calculation
masknoise=1;
% threshold for noise
thold=3; % signal has to be thold * mean(noise) to be used!!

% Filter the data before doing anything further:
DATA{1}.FTDATA=l_FTfilter(DATA{1}.FTDATA);


bg_roi=[];
if ~isempty(ROI)
    for ii=1:length(ROI)
        % assume background roi is labeled 'bg'
        if strcmp(ROI(ii).label,'bg')
            bg_roi=ii;
        end
    end

    
    % now there either is bg-roi or not..
    % separate slices:
    if ns>1
        for ii=1:ns;
            d(ii).data=double(DATA{1}.FTDATA(:,:,ii,:));  % slices in #3
            d(ii).T1=double(T1DATA{1}.FTDATA(:,:,ii));  % slices in #3
        end
    else
        d(1).data=DATA{1}.FTDATA;
    end
    if ~isempty(bg_roi)
        % also bg-ROIs present.. grab noise-estimate vector:
        noise=double(DATA{1}.FTDATA(ROI(bg_roi).voxels{1}(:,:,1)));
        % save noiselevel:
        NoiseLevel=mean(noise);

    end
end

if isempty(bg_roi)
    % no noise-ROI found, set noise level to zero (ie, do thresholding)
    noise=0;
end


%keyboard;

%% check if name was given, then don't ask
if (nargin<5) || isempty(savefname)
% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.t1r;*.T1R;*.s1r;*.S1R;*.MAT;*.mat',...
                    'T1R-Files (*.t1r, *.s1r, *.mat)';...
                    '*.*','All Files (*.*)'},...
                               'Save T1R-file',[DATA{1}.HDR.fpath, ...
                    'SE0000_t1rho_map.mat']);
if isequal(fname,0) || isequal(fpath,0)
    f=[];
    return
end
else
    % save-filename was given, make compatible with the rest of the
    % function
    [fpath,fname,fext]=fileparts(savefname);
    fname=[fname fext];
end



%keyboard;
fitopts=optimset('fminsearch');
fitopts.Display='off';

% FIND CALCULATION ROI
calc_roi=[];
if ~isempty(ROI)
    for ii=1:length(ROI)
        if strcmp(ROI(ii).label,'ROI 1') ||strcmp(ROI(ii).label,'calc_roi') ||strcmp(ROI(ii).label,'calc roi') 
            calc_roi=ii;
        end
    end
end

%keyboard;
tic


% Loop through slices:
for nn=1:ns
 % keyboard;  
    % USE ONLY ROI AREA!!
    try
    if ~isempty(calc_roi)
        % assume this ROI is found at the same slice as the image...:
        RR=ROI(calc_roi).voxels{1}(:,:,nn);
        
        if ns>1
            % if multislice, repmat by #4
            RR=repmat(RR,[1 1 1 size(d(nn).data,4)]);
        else
            % if single, repmat by #3 
            RR=repmat(RR,[1 1 size(d(nn).data,3)]);
        end
        
        d(nn).data=squeeze(d(nn).data.*RR);  % If existing, remove the empty 3rd dim.

        if masknoise
            % Find image with highest signal - and compare to noise
            tmp=d(nn).data;
            tmp=reshape(tmp,[size(tmp,1)*size(tmp,2),size(tmp,3)]); % last dim has different TE/TR/whatev
            [mval,mind]=max(mean(tmp));
            
            % use mean as noise value:
            noise_limit=mean(noise)*thold;
            % or use std as noise value:
            %noise_limit=std(noise)*thold;
            
            % save noiselevel:
            NoiseLevel=mean(noise);

        
            % create mask and use it:
            mask=ones(size(d(nn).data(:,:,1)));
            % zero all pixels where noise_limit is hit:
            mask(d(nn).data(:,:,mind)<noise_limit)=0;

            % grow mask to data size
            mask=repmat(mask,[1,1,size(d(nn).data,3)]);
            d(nn).data=d(nn).data.*mask;
        end

    else
        % do nothing..
    end
    catch
     %   keyboard;
        errordlg({sprintf('Error in masking noise (step %1.0d). The following error was returned',ii),...
            '',lasterr},'modal')
    end



    % Calculate the map
    [fp,fn,fe]=fileparts([fpath,fname]);
    

   % try
        m=size(d(nn).data,1);n=size(d(nn).data,2);
        modder=floor(m/50);
        S0=[];
        RR=[];
        Err=[];
        Baseline=[];
    
        h=waitbar(0,sprintf('calculating slice %1.0d / %1.0d',nn,ns));
        drawnow;

        % try to speed up calculation by permuting data beforehand for
        % easier matrix notation:
        tmp=double(permute(d(nn).data,[3 1 2]));
        if length(nex)==size(tmp,3)
            % if variable number of averages used, divide them out:
            nexmat=repmat(permute(nex(:),[2 3 1]),[size(tmp,1) size(tmp,2) 1]);
            tmp=tmp./nexmat;
        end
        
        outputsize=[3 3 4 4];
        for ii=1:m
            for jj=1:n
                s=tmp(:,ii,jj);

                
                if (~sum(s)==0)
                    if fit_type==4
                        % fit with T1
                        t1=d(nn).T1(ii,jj);
                        %keyboard;

                        
                        th2=l_calculate_t2(tSL,s,fit_type,0,NoiseLevel,t1);
                        
                        % ------TODO--------
                        %  T1-fitting! All the data is here:
                        %  tSL, s, T1
                        
                    else
                        th2=l_calculate_t2(tSL,s,fit_type,0,NoiseLevel);
                    end                    
                else
                    th2=zeros(1,outputsize(fit_type));
                end
                
                % save values
                S0(ii,jj)=th2(1);
                % T1r-map
                RR(ii,jj)=th2(2);
  
                if outputsize(fit_type)==3
                    % fitting error
                    Err(ii,jj)=th2(3);
                elseif outputsize(fit_type)==4
                    Baseline(ii,jj)=th2(3); % baseline
                    % fitting error
                    Err(ii,jj)=th2(4);
                end
            end
            if mod(ii,modder)==0
                waitbar(ii/m,h);
            end
        end
        close(h);
    
  
        % Save data
        Param(nn).S0=S0;
        Param(nn).FitValues=tSL;
        Param(nn).Type='Uncorrected Adiabatic T1rho';
        Param(nn).Linear=0;
        Param(nn).Err=Err;
        Param(nn).RR=RR;
        Param(nn).BL=Baseline;
        Param(nn).Range=[0 3000];
        
        % LIMIT DATA RANGE!!!!!!
        RR(RR>3000)=0;
        RR(RR<0)=0;
        
        Data(:,:,nn)=RR;
        

if isfield(DATA{1}.HDR,'dicominfo')
    % DICOM FOUND, USE IT!
        % stuff the dicom headers in for later use:
        Param(nn).Header=DATA{1}.HDR.dicominfo(nn:ns:end);
end

    
   % catch
   %     errordlg({'Could not calculate maps. The following error was returned',...
   %         '',lasterr},'modal')
   % end

end

clear DATA;

% Generate a "More Aedes Compatible" file:
DATA.FTDATA=Data;
DATA.Param=Param;
DATA.DataFormat='mat';
DATA.HDR.fname=fname;
DATA.HDR.fpath=[fpath filesep];

if exist('p')
    HDR=p;
    DATA.HDR.dicominfo=p;
end

f=DATA;

% All slices calculated, save data
    save(fullfile(fpath,fname),'DATA');
toc





function f=l_fit_MT(DATA,ROI,fit_type,T1DATA,savefname,resp)
% This Aedes plugin calculates MT/T1sat map from Siemens-MT-prep data

% This function is written for Aedes
%
% copyright (c) 2013 Mikko Nissi <nissi@cmrr.umn.edu>
%
% Changes
% 07/2013, Mikko Nissi
%    - Initial version
%    - Supports 4D-data containing both inversion-prepared and non-prepared
%    volumes
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


% check input for T1DATA
if (nargin == 6)
    clearvars T1DATA;
    T1DATA{1}.FTDATA=zeros(size(DATA{1}.FTDATA(:,:,:,1)));
end
if nargin<4
    T1DATA{1}.FTDATA=zeros(size(DATA{1}.FTDATA(:,:,:,1)));   
end


n_inv=1;  % default dummy value.

if isfield(DATA{1}.HDR,'dicominfo')
    % DICOM FOUND, USE IT!

    
    % Grab dicom stuff
    p=DATA{1}.HDR.dicominfo;
    
    % Dig TR from PROCPAR (in seconds)
    TR = p(1).RepetitionTime;
    
    % Check size of dicomheader and guess the data properties based on this
    if size(DATA{1}.HDR.dicominfo,1)==2
        % two sets of dicomfiles in the data
        two_vol=1;
    else
        two_vol=0;
    end
    
    N_volumes=size(DATA{1}.HDR.dicominfo,1);
    
    % Check which volume is inverted if multiple
    if N_volumes==2
        [fl1,fd1,strrs1]=l_new_parse_swip(DATA,1); % first volume checked
        [fl2,fd2,strrs2]=l_new_parse_swip(DATA,2); % second volume checked
        if fl1(5)==fl2(5)
            disp('Both volumes with or without inversion. Not supported');
            f=[];
            return;
        end
        % ok, volumes differ, one has to be inverted
        [void,n_inv]=min([fl1(5) fl2(5)]); % fl(5)==1 => -Z, fl(5)==2 => +Z
    else
        n_inv=1; % default dummy value to make this work if only one volume used
    end
    
    
    % dig sWiPMemBlock of the first header
    [fl,fd,strrs]=l_new_parse_swip(DATA,1);
    
    % WIP-data, naming follows sequence internal variable names
    WIP_RectPulseDur = fl(1);
    WIP_NumPrepPulses= fl(2);
    WIP_PrepPulseDur= fl(3);
    WIP_T2rhoPulseDur= fl(4);
    WIP_InversionCheckBox= fl(5);
    WIP_InvPulseDur= fl(6);
    WIP_SpoilerMom= fl(7);
    WIP_MTOffsetFreq= fl(8);
    WIP_T2rhoCheckBox= fl(9);
    WIP_ArrayCheckBox= fl(10);
    WIP_InternalTR= fl(11);
    WIP_PrepPulseShape= fl(12);   % CHECK THIS FOR RAFF or MT!
    WIP_RAFFPulseDur= fl(13);
    WIP_NumRectPulses= fl(14);
    %WIP_MinIPD= fl(15);
    
    dWIP_PrepVoltageAdj = fd(1);
    dWIP_RAFFVoltageAdj= fd(2);
    dWIP_RectVoltageAdj= fd(3);
    dWIP_T2rVoltageAdj= fd(4);
    dWIP_InvVoltageAd= fd(5);
    
    
    % Number of prep pulses per block
    % Hard-coded in the sequence at the moment!!! FIXTHIS!!
    NPPB=4;
    
    
    if two_vol
        % dig pars from the second header (assumed to be from the second data
        % set because it will be loaded with d_loader. Or otherwise this
        % function has no meaning anyways..)
        [fl2,fd2,strrs2]=l_new_parse_swip(DATA,2);
        rp=[2 5 7 10 12 13];   % Relevant WIP parameters to check!
        
        % do some checks..
        if sum(abs(fl(rp)-fl2(rp)))==0 || sum(abs(fl(rp)-fl2(rp)))>1
            % something is not right, other than inversion pulse changed. abort
            disp('Measurements do not match. Check your data acquisition.');
            f=[];
            return;
        end
        % this should give mask for "isinverted", because '1' indicates yes and
        % '2' indicates 'no':
        aa=abs([fl(5) fl2(5)]-2);
        n_inv = find(aa); % index of the inverted data
        % presumably data is in order
    end
    
    
    if WIP_PrepPulseShape==3 || WIP_PrepPulseShape==4|| WIP_PrepPulseShape==5|| WIP_PrepPulseShape==6|| WIP_PrepPulseShape==7|| WIP_PrepPulseShape==8
        % It's a RAFF SEQUENCE!
        disp(sprintf('Prep shape = %1.0f, fitting for RAFF',WIP_PrepPulseShape));
        pw=WIP_RAFFPulseDur;  % us, RAFF pulse duration
        pw=pw/1000; % ms
        
        % For now, the array of pulses is as follows:
        if (WIP_ArrayCheckBox==1)
            n_mlev=[0:NPPB:WIP_NumPrepPulses];
        end
    else
        % It's an MT SEQUENCE!
        disp(sprintf('Prep shape = %1.0f, fitting for T1SAT (MT)',WIP_PrepPulseShape));
        pw=WIP_RectPulseDur;  % RECT pulse duration, in ms
        
        % For now, the array of pulses is as follows:
        if (WIP_ArrayCheckBox==1)
            n_mlev=[0:NPPB:WIP_NumRectPulses];
        end
    end
    
    
    % create tSL-vector
    tSL=pw.*n_mlev
    
    
    
else
    % NO DICOM.. display dialog
    if(nargin == 3)
    resp = aedes_inputdlg('Type tSL values (just once even if +Z & -Z data)','Input dialog',...
        {mat2str([])});
    end
    if isempty(resp)
        % NO RESPONSE OR SOME OTHER ISSUES.. EXIT
        f=[];
        return
    else
        resp=resp{1};
        fit_vals = str2num(resp);
        tSL=fit_vals;
        if size(DATA{1}.FTDATA,ndims(DATA{1}.FTDATA))>length(tSL)
            % presume two volume dataset
            n_inv=2;
        end
        
    end
end


%keyboard;

% averages: this has no meaning here..
nex=1;

% number of slices - from d_loader, the volumes contain TE or TR /whatnot
% and the dim(3) contains slices, IF there are 4 dimensions to the data
if (ndims(DATA{1}.FTDATA))==4
    if size(DATA{1}.FTDATA,4)==2
        % this is most likely single-slice data
        ns=1;
        
        % RESHAPE DATA TO THE SAME FORMAT AS MULTISLICE:
        % [x, y, slices, tSL]
        szm=size(DATA{1}.FTDATA);
        DATA{1}.FTDATA=reshape(DATA{1}.FTDATA,[szm(1),szm(2),1,szm(3)*szm(4)]);
    end
    % Now dim3 should contain correct number of slices
    ns=size(DATA{1}.FTDATA,3);
else
    ns=1;
end



%keyboard; %% find background ROI
% if masknoise, then signals below noise-threshold are not considered for
% calculation
masknoise=1;
% threshold for noise
thold=3; % signal has to be thold * mean(noise) to be used!!

% Filter the data before doing anything further:
DATA{1}.FTDATA=l_FTfilter(DATA{1}.FTDATA);
%keyboard;

bg_roi=[];
if ~isempty(ROI)
    for ii=1:length(ROI)
        % assume background roi is labeled 'bg'
        if strcmp(ROI(ii).label,'bg')
            bg_roi=ii;
        end
    end
    
    
    % now there either is bg-roi or not..
    % separate slices:
    if ns>1
        for ii=1:ns;
            d(ii).data=double(DATA{1}.FTDATA(:,:,ii,:));  % slices in #3
            d(ii).T1=double(T1DATA{1}.FTDATA(:,:,ii));  % slices in #3
        end
    else
        d(1).data=DATA{1}.FTDATA;
    end
    if ~isempty(bg_roi)
        % also bg-ROIs present.. grab noise-estimate vector:
        noise=double(DATA{1}.FTDATA(ROI(bg_roi).voxels{1}(:,:,1)));
        % save noiselevel:
        NoiseLevel=mean(noise);
        
    end
else
    % No ROIs whatsoever..
    % separate slices:
    if ns>1
        for ii=1:ns;
            d(ii).data=double(DATA{1}.FTDATA(:,:,ii,:));  % slices in #3
            d(ii).T1=double(T1DATA{1}.FTDATA(:,:,ii));  % slices in #3
        end
    else
        d(1).data=DATA{1}.FTDATA;
    end
    NoiseLevel=0;

end

if isempty(bg_roi)
    % no noise-ROI found, set noise level to zero (ie, do thresholding)
    noise=0;
end


%keyboard;

%% check if name was given, then don't ask
if (nargin<5) || isempty(savefname)
% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.t1sat;*.T1SAT;*.s1sat;*.S1SAT;*.MAT;*.mat',...
                    'T1SAT-Files (*.t1sat, *.s1sat, *.mat)';...
                    '*.*','All Files (*.*)'},...
                               'Save T1SAT-file',[DATA{1}.HDR.fpath, ...
                    't1sat_map.mat']);
if isequal(fname,0) || isequal(fpath,0)
    f=[];
    return
end
else
    % save-filename was given, make compatible with the rest of the
    % function
    [fpath,fname,fext]=fileparts(savefname);
    fname=[fname fext];
end



%keyboard;
fitopts=optimset('fminsearch');
fitopts.Display='off';

% FIND CALCULATION ROI
calc_roi=[];
if ~isempty(ROI)
    for ii=1:length(ROI)
        if strcmp(ROI(ii).label,'ROI 1') ||strcmp(ROI(ii).label,'calc_roi') ||strcmp(ROI(ii).label,'calc roi') 
            calc_roi=ii;
        end
    end
end

%keyboard;
tic


% Loop through slices:
for nn=1:ns
 % keyboard;  
    % USE ONLY ROI AREA!!
    try
    if ~isempty(calc_roi)
        % assume this ROI is found at the same slice as the image...:
        RR=ROI(calc_roi).voxels{1}(:,:,nn);
        
        %if ns>1
            % if multislice, repmat by #4
            RR=repmat(RR,[1 1 1 size(d(nn).data,4)]);
        %else
            % if single, repmat by #3 
        %    RR=repmat(RR,[1 1 size(d(nn).data,3) size(d(nn).data,4)]);
        %end
        
        d(nn).data=squeeze(d(nn).data.*RR);  % If existing, remove the empty 3rd dim.

        if masknoise
            % Find image with highest signal - and compare to noise
            tmp=d(nn).data;
            tmp=reshape(tmp,[size(tmp,1)*size(tmp,2),size(tmp,3)]); % last dim has different TE/TR/whatev
            [mval,mind]=max(mean(tmp));
            
            % use mean as noise value:
            noise_limit=mean(noise)*thold;
            % or use std as noise value:
            %noise_limit=std(noise)*thold;
            
            % save noiselevel:
            NoiseLevel=mean(noise);

        
            % create mask and use it:
            mask=ones(size(d(nn).data(:,:,1)));
            % zero all pixels where noise_limit is hit:
            mask(d(nn).data(:,:,mind)<noise_limit)=0;

            % grow mask to data size
            mask=repmat(mask,[1,1,size(d(nn).data,3)]);
            d(nn).data=d(nn).data.*mask;
        end

    else
        % do nothing..
        % Not actually nothing -- still remove the extra dim!
        d(nn).data=squeeze(d(nn).data);  % If existing, remove the empty 3rd dim.
                
    end
    catch
     %   keyboard;
        errordlg({sprintf('Error in masking noise (step %1.0d). The following error was returned',nn),...
            '',lasterr},'modal')
    end



    % Calculate the map
    [fp,fn,fe]=fileparts([fpath,fname]);
    
%keyboard;

   % try
        m=size(d(nn).data,1);n=size(d(nn).data,2);
        modder=floor(m/50);
        S0=[];
        RR=[];
        SSS=[];
        Sloss=[];
        Resid=[];
    
        h=waitbar(0,sprintf('calculating slice %1.0d / %1.0d',nn,ns));
        drawnow;

        % try to speed up calculation by permuting data beforehand for
        % easier matrix notation:
        tmp=double(permute(d(nn).data,[3 1 2]));
        if length(nex)==size(tmp,3)
            % if variable number of averages used, divide them out:
            nexmat=repmat(permute(nex(:),[2 3 1]),[size(tmp,1) size(tmp,2) 1]);
            tmp=tmp./nexmat;
        end

        

                
        
        outputsize=[4 5 4 5];
        fitchoose=[3 4 5 6];
        
        for ii=1:m
            for jj=1:n
                s=tmp(:,ii,jj);

                if (~sum(s)==0)
                    if fit_type==4 || fit_type==3
                        % fit with T1
                        t1=d(nn).T1(ii,jj);
                        th2=l_calculate_traff(tSL,s,0,fitchoose(fit_type),n_inv,NoiseLevel,t1);
                        
                        % ------TODO--------
                        %  T1-fitting! All the data is here:
                        %  tSL, s, T1
                        
                    else                    
                        % fit RAFF/MT using given number of parameters..
                        th2=l_calculate_traff(tSL,s,0,fitchoose(fit_type),n_inv,NoiseLevel); 
                    end
                else
                    th2=zeros(1,outputsize(fit_type));
                end
            
                % Grab fitting values..
                S0(ii,jj)=th2(1);  % S_0
                RR(ii,jj)=th2(2);  % Traff
                SSS(ii,jj)=th2(3); % S_ss

                if outputsize(fit_type)==4
                    % fitting error
                    Resid(ii,jj)=th2(4); % fitting residual
                elseif outputsize(fit_type)==5
                    % If 4-parameter fit:
                    Sloss(ii,jj)=th2(4); % "inversion loss"
                    % fitting error
                    Resid(ii,jj)=th2(5); % fitting residual
                end
                
            end
            if mod(ii,modder)==0
                waitbar(ii/m,h);
            end
        end
        close(h);
    
        % Save data
        Param(nn).S0=S0;
        Param(nn).FitValues=tSL;
        Param(nn).Type='T1sat';
        Param(nn).Linear=0;
        Param(nn).Sss=SSS;
        Param(nn).RR=RR;
        Param(nn).Sloss=Sloss;  % may be empty..
        Param(nn).Err=Resid;
        Param(nn).Range=[0 3000]; % dummy estimate to aid later dicom export..
        
        % LIMIT DATA RANGE!!!!!!
        RR(RR>3000)=0;
        RR(RR<0)=0;

        Data(:,:,nn)=RR;

      

if isfield(DATA{1}.HDR,'dicominfo')
    % DICOM FOUND, USE IT!
        % stuff the dicom headers in for later use:
        Param(nn).Header=DATA{1}.HDR.dicominfo(nn:ns:end);
end

   
    
   % catch
   %     errordlg({'Could not calculate maps. The following error was returned',...
   %         '',lasterr},'modal')
   % end

end

clear DATA;

% Generate a "More Aedes Compatible" file:
DATA.FTDATA=Data;
DATA.Param=Param;
DATA.DataFormat='mat';
DATA.HDR.fname=fname;
DATA.HDR.fpath=[fpath filesep];


if exist('p')
    HDR=p;
    DATA.HDR.dicominfo=p;
end


f=DATA;

% All slices calculated, save data
    save(fullfile(fpath,fname),'DATA');

toc




function [fl,fd,strings]=l_new_parse_swip(DATA,num)

if nargin<2
    num=1;
end

hdr=DATA{1}.HDR;
dcm=hdr.dicominfo(num);
stuff=char(dcm.Private_0029_1020)';
stuff2=stuff(strfind(stuff,'### ASCCONV BEGIN'): strfind(stuff,'### ASCCONV END'));
rr=textscan(stuff2,'%s\n',10000,'Delimiter','');


%keyboard;
gerp=rr{1};
per={};
for ii=1:length(gerp)
    if regexpi(upper(gerp{ii}),'SWIPMEM')
        per{end+1}=gerp{ii};
    end
end

% grep doesn't exist in matlab commands..
%per=grep(rr{1},'-x','sWiPMem');

% keyboard
fl=[];
fd=[];
for ii=1:length(per)
    tmpvals=sscanf(upper(per{ii}),'SWIPMEMBLOCK.ALFREE[%f] = %f');
    if ~isempty(tmpvals)
        fl(end+1,:)=tmpvals(:)';
    end
    tmpvads=sscanf(upper(per{ii}),'SWIPMEMBLOCK.ADFREE[%f] = %f');
    if ~isempty(tmpvads)
        fd(end+1,:)=tmpvads(:)';
    end
end

inc1=1-fl(1,1);
inc2=1-fd(1,1);

tmp1(inc1+round(fl(:,1)))=fl(:,2);
tmp2(inc2+round(fd(:,1)))=fd(:,2);


%[val,sor]=sort(fl(:,1));
%fl=fl(sor,2);
%[val,sor]=sort(fd(:,1));
%fd=fd(sor,2);
fl=tmp1;
fd=tmp2;
strings=per(:);





function f=l_FTfilter(imdata,fstrength,sw)
%Fourier (low-pass) filter 2D image or image stack
%    X = FTFILTER(IMDATA,FSTRENGTH,SWITCH) filters image data in IMDATA
%    using filtering strength FSTRENGTH. SWITCH, if set, does an additional
%    histogram filtering (may help if IMDATA has spikes in it, otherwise
%    better not use).
%

% (c) Mikko Nissi, 2012
% Based on original code by Timo Liimatainen

if nargin<3
    sw=0;
end
if nargin<2
    fstrength=.3;
end

%filtering data
filterweight=hanning(size(imdata(:,:,1),1))*hanning(size(imdata(:,:,1),2))';
filterweight=filterweight.^fstrength;
    
for ii=1:length(imdata(1,1,:))
    if sw
    % First do histogram-equalization-filtering
        imdata(:,:,ii)=histofilter(imdata(:,:,ii),0.01,[]); % remove data until histogram min is within 0.01 of max
    end
    % Filter data
    imdata(:,:,ii)=abs(fft2(ifftshift(fftshift(ifft2(imdata(:,:,ii))).*filterweight)));
end

f=imdata;


%
%
%
%
%
%
%
%   FITTING FUNCTIONS 
%
%
%
%
%
%
%
%
%

function f=l_calculate_t2(TE,signal,npars,plotfit,noiselevel,t1)
%CALCULATE_T2 Calculates T2 for given TE and signal
%    CALCULATE_T2(TE,SIGNAL,NPARS)
%    Outputs M0, T2 and fitting error with NPARS=2 and
%    M0, T2, Constant and fitting error with NPARS=3

% 8/2013, Mikko Nissi, CMRR
%  - Changed error to be scaled by the maximum of the signal to provide
%    more meaningful estimate of the goodness of the fit! Probably
%    reasonable to estimate, that if error (2-norm of diff) < 0.5 x max,
%    the fit is ok, otherwise not.
%  - Added similar norm-diff as the fitting error of linearized fit
%
% 3/2012, Mikko Nissi, CMRR
%  - Added noise level consideration to fitting
%
% 31.8.2010, Mikko Nissi, UEF
%  - Fixed figure opening with plotfit=0..
% 30.8.2010, Mikko Nissi, UEF
%  - Modified fitting functions -> S0 separated for different T2s as
%    it should be..
% 20.4.2010, Mikko Nissi, Uni. E. Finland
%  - comments added..
%  - plotting by parameter.. may slow down due to one more check on every
%  run
% 17.11.2008, revised fitting functions..
% (c) Mikko Nissi 4.11.2008, Kuopio university hospital and University of
% Kuopio
%




%% Check input
if(nargin<2)
    disp('Need data to work...');
    return;
end;

if nargin<3
    % default to 1 peak + baseline
    npars=3;
end;

if nargin<4
    % default to no plotting..
    plotfit=0;
else
    if plotfit==1
        hh=figure;
    end
end

if nargin<5
    % default to noise equal 0
    noiselevel=0;
elseif isempty(noiselevel)
    noiselevel=0;
end


TE=TE(:);
signal=signal(:);

% fminsearch options
fitpars=optimset('fminsearch');
%fitpars.MaxFunEvals=100;
%fitpars.MaxIter=100;
%fitpars.TolFun=1;
fitpars.Display='off';

%% Check fit type..
switch npars
    case 1
        % linearize!
        H=[-TE ones(size(TE))];
        if (min(signal)~=0)
            s=log(signal);
            th2=pinv(H'*H)*H'*s;
        else
            th2=[1 1];
        end;
        th2=[exp(th2(2)) 1/th2(1)];
        t2est=th2(1)*exp(-TE/th2(2));
        erhor=norm(t2est(:)-signal(:),2);
        
    case 2 % two parameters: s0, T2
        th0=[max(signal) max(TE)];
        [th2,erhor]=fminsearch(@t2_2par,th0,fitpars,TE,signal,noiselevel,plotfit);
    case 3 % three parameters: s0, T2, baseline
        th0=[max(signal) max(TE) noiselevel];
        [th2,erhor]=fminsearch(@t2_3par,th0,fitpars,TE,signal,noiselevel,plotfit); 
    case 4 % four parameters: two peaks, no baseline
        th0=[max(signal) min(TE) mean(signal) mean(TE)];
        [th2,erhor]=fminsearch(@t2_2peak,th0,fitpars,TE,signal,noiselevel,plotfit); 
    case 5 % five parameters: two peaks + baseline
        th0=[max(signal) min(TE) mean(signal) mean(TE) noiselevel];
        [th2,erhor]=fminsearch(@t2_3par_2peak,th0,fitpars,TE,signal,noiselevel,plotfit); 
end

% Scale error with the signal maximum:
erhor=erhor/max(signal);

% Try empirically zero ill-fitting (bad) data:
if erhor>1
    th2=zeros(size(th2));
end

% return values:        
f=[th2 erhor];
return;


%% Fitting functions

function f=t2_2par(th0,te,orig,noiselevel,plotfit) % case 2: 1 peak, no baseline

  y=th0(1)*exp(-te./th0(2))+noiselevel;
  f=norm(y(:)-orig(:),2);
  
  if plotfit
      % Constants
      linewid=1.5;
      markersz=10;
      tt=linspace(0,max(te)*1.1,50);
      yy=th0(1)*exp(-tt./th0(2))+noiselevel;
      plot(te,orig,'bx',tt,yy,'r:','linewidth',linewid,'markersize',markersz);drawnow;
  end
  return;

function f=t2_3par(th0,te,orig,noiselevel,plotfit) % case 3: 1 peak, baseline
  y=th0(1)*exp(-te./th0(2))+th0(3);
  f=norm(y(:)-orig(:),2);
  if plotfit
      % Constants
      linewid=1.5;
      markersz=10;

      tt=linspace(0,max(te)*1.1,50);
      yy=th0(1)*exp(-tt./th0(2))+th0(3);
      plot(te,orig,'bx',tt,yy,'r:','linewidth',linewid,'markersize',markersz);drawnow;
  end
  return;

function f=t2_2peak(th0,te,orig,noiselevel,plotfit) % case 4: 2 peaks, no baseline
  y=th0(1)*exp(-te./th0(2))+th0(3)*exp(-te./th0(4))+noiselevel;
  f=norm(y(:)-orig(:),2);
  if plotfit
      % Constants
      linewid=1.5;
      markersz=10;

      tt=linspace(0,max(te)*1.1,50);
      yy=th0(1)*exp(-tt./th0(2))+th0(3)*exp(-tt./th0(4))+noiselevel;
      plot(te,orig,'bx',tt,yy,'r:','linewidth',linewid,'markersize',markersz);drawnow;
  end

function f=t2_3par_2peak(th0,te,orig,noiselevel,plotfit) % case 5: 2 peaks, baseline
  y=th0(1)*exp(-te./th0(2))+th0(3)*exp(-te./th0(4))+th0(5);
  f=norm(y(:)-orig(:),2);
  if plotfit
      % Constants
      linewid=1.5;
      markersz=10;

      tt=linspace(0,max(te)*1.1,50);
      yy=th0(1)*exp(-tt./th0(2))+th0(3)*exp(-tt./th0(4))+th0(5);
      plot(te,orig,'bx',tt,yy,'r:','linewidth',linewid,'markersize',markersz);drawnow;
  end









function f=l_calculate_traff(tSL,S,plotfit,tyyp,n_inv,NoiseLevel,t1)
%CALCULATE_TRAFF Calculates T_RAFF relaxation time constant using input values
%    TH=CALCULATE_TRAFF(TSL,S,PLOTF,TYPE,N_INV)
%    Returns: S_0, T_RAFF, S_SS [,inversion loss] and FIT_RESIDUAL
%    - TSL contains pulse durations for each preparation
%    - S contains measured signal
%    - PLOTF is a flag whether or not to plot the curve fitting
%    - TYPE is  the number of parameters (3 or 4)
%    - N_INV indicates which half of the signal was inversion-prepared (1 or 2)

% Revisions
% 2011-2013, Mikko Nissi, CMRR
%    - Added 4-parameter fitting
%    - Some bug fixes and improved guesswork
%
%

% (c) Mikko Nissi 2009, University of Kuopio, based on code by
% Timo Liimatainen. See Liimatainen et al, MRM 2010


%% Check input
if(nargin<2)
    disp('Need data to work...');
    return;
end;

if nargin<3
    % default to no plotting..
    plotfit=0;
else
    if plotfit==1
        hh=figure;
    end
end

if nargin<4
    % default to 3-parameter fit
    tyyp=3;
end

if nargin<5
    % number of inverted volume not given, assume SECOND
    n_inv=2;
end
    
    
tSL=tSL(:);
S=S(:);

% fminsearch options
fitpars=optimset('fminsearch');
%fitpars.MaxFunEvals=100;
%fitpars.MaxIter=100;
%fitpars.TolFun=1;
fitpars.Display='off';


if tyyp==4
    % Initial guess
    th0=[max(S) max(tSL)/2 min(S) max(S)/10];
else
    % Initial guess
    th0=[max(S) max(tSL)/2 min(S)]; % max(S)/10];
end

% Actually inversion can be guessed if length(TSL)==length(S)/2.. Use this!
if length(tSL)*2==length(S);
    % inverted, invert ONE half of the data, either indicated or second half:
    S((end/2)*(n_inv-1)+1:(end/2)*(n_inv))= -S((end/2)*(n_inv-1)+1:(end/2)*(n_inv));
    if tyyp==4
        [th2,erhor]=fminsearch(@traff4_inv,th0,fitpars,tSL,S,plotfit);
    else
        [th2,erhor]=fminsearch(@traff_inv,th0,fitpars,tSL,S,plotfit);
    end        
else
    % not inverted
    [th2,erhor]=fminsearch(@traff,th0,fitpars,tSL,S,plotfit);
end
        
% return values:
f=[th2 erhor];
return;




function f=traff_inv(th0,TSL,S,plotfit)

if nargin<4
    plotfit=0;
end
% th0(1) = S0
% th0(2) = T_RAFF
% th0(3) = SS
% th0(4) = inversion loss..

y=[th0(1)*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2))); ...
   -(th0(1))*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2)))];
%   -(th0(1)-th0(4))*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2)))];


f=norm((y(:))-S(:),2);

if plotfit==1
    tt=linspace(0,max(TSL)*1.1,50);
    yy1=th0(1)*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
    yy2=-(th0(1))*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
  %     yy2=-(th0(1)-th0(4))*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
    
    plot(TSL,S(1:floor(end/2)),'bo',TSL,S(floor(end/2)+1:end),'ro',tt,yy1,'g:',tt,yy2,'k:','markersize',8,'linewidth',2);
    drawnow;
end

return;




function f=traff4_inv(th0,TSL,S,plotfit)

if nargin<4
    plotfit=0;
end
% th0(1) = S0
% th0(2) = T_RAFF
% th0(3) = SS
% th0(4) = inversion loss..

y=[th0(1)*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2))); ...
   -(th0(1)-th0(4))*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2)))];
%   -(th0(1))*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2)))];


f=norm((y(:))-S(:),2);

if plotfit==1
    tt=linspace(0,max(TSL)*1.1,50);
    yy1=th0(1)*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
  %  yy2=-(th0(1))*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
       yy2=-(th0(1)-th0(4))*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
    
    plot(TSL,S(1:floor(end/2)),'bo',TSL,S(floor(end/2)+1:end),'ro',tt,yy1,'g:',tt,yy2,'k:','markersize',8,'linewidth',2);
    drawnow;
end

return;






%% Non-inverted version
function f=traff(th0,TSL,S,plotfit)

if nargin<4
    plotfit=0;
end
% th0(1) = S0
% th0(2) = T_RAFF
% th0(3) = SS

y=th0(1)*exp(-TSL./th0(2))+th0(3)*(1-exp(-TSL./th0(2)));


f=norm((y(:))-S(:),2);

if plotfit==1
    tt=linspace(0,max(TSL)*1.1,50);
    yy1=th0(1)*exp(-tt./th0(2))+th0(3)*(1-exp(-tt./th0(2)));
    plot(TSL,S,'bo',tt,yy1,'g:','markersize',8,'linewidth',2);
    drawnow;
end

return;

function varargout = inverse_cat(DIM,C)
% INVERSE_CAT splits data into sub-arrays along the specified dimension.
%     [A B]=INVERSE_CAT(DIM,M) splits array M along dimension, DIM,
%     returning sub-arrays A and B.
%
%  Examples:
%     M = [1 2 3; 4 5 6; 7 8 9];
%     C = cat(2,M,M)
%     [A B] = inverse_cat(2,C)         ... returns A=M and B=M
%     [A B] = inverse_cat(1,rot90(C))  ... returns A=rot90(M) B=rot90(M)
%     [A B] = inverse_cat(3,cat(3,M,M))... returns A=M and B=M
%
%     See also cat, num2cell.

if (~(DIM==round(DIM)) || DIM<1)
    error('Dimension must be a finite integer.');
end

nout = max(nargout,1);
x = cell(nout,1);
q = floor(size(C,DIM)/nout);

if (isempty(C)), varargout = x; return; end

%separate along dimension 2
if (DIM==1)
    for i=1:nout
    	x{i} = C(i*q - q + 1:i*q, :);
    end
    %include extra data from uneven division
    if (mod(size(C,DIM),nout)~=0)
        x{nout} = C(nout*q - q + 1:end, :);
    end
%separate along dimension 2
elseif (DIM==2)
    for i=1:nout
        x{i} = C(:, i*q - q + 1:i*q);
    end
    %include extra data from uneven division
    if (mod(size(C,DIM),nout)~=0)
        x{nout} = C(:, nout*q - q + 1:end);
    end
%separate along dimension 3
elseif (DIM==3)
    for i=1:nout
        x{i} = C(: , :, i*q - q + 1:i*q);
    end
    %include extra data from uneven division
    if (mod(size(C,DIM),nout)~=0)
        x{nout} = C(:, :, nout*q - q + 1:end);
    end
%separate along arbitrary dimension
else 
    index = cell(1, ndims(C)); 
    index(:) = {':'}; 
    for i = 1:nout 
      index{DIM} = i*q - q + 1:i*q; 
      x{i} = C(index{:}); 
    end
end

varargout=x;



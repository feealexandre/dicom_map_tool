function export_dicom(DATA,ROI,AddInfo)
% This Aedes plugin exports FTDATA as dicom file

% This function is written for Aedes
%
% Updates
% 2/26/2016, Mikko Nissi
%    - Additional check for parametric data if not recognized origin
%
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
         %   keyboard;
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


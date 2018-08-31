function [fl,fd,strings]=new_parse_swip(DATA,num)
%Parses Siemens DICOM WIP memory block to vectors of values [and strings]

% Mikko Nissi, 2012

if nargin<2
    num=1;
end

hdr=DATA{1}.HDR;

% dcm to contain the DICOM HEADER of the image of interest, currently
% grabbed from Aedes-DATA structure.
dcm=hdr.dicominfo(num);
stuff=char(dcm.Private_0029_1020)';
stuff2=stuff(strfind(stuff,'### ASCCONV BEGIN'): strfind(stuff,'### ASCCONV END'));
rr=textscan(stuff2,'%s\n',10000,'Delimiter','');

gerp=rr{1};
per={};
for ii=1:length(gerp)
    if regexpi(upper(gerp{ii}),'SWIPMEM')
        per{end+1}=gerp{ii};
    end
end

% Note: grep doesn't exist in matlab commands..
%per=grep(rr{1},'-x','sWiPMem');

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

% Use the vector indices from SWIPMEMBLOCK for output vectors (keep the
% same WIP numbering even if some of these are missing/missed for some
% reason in the above sscanf)
inc1=1-fl(1,1);   % offset to 1, just in case
inc2=1-fd(1,1);   % 
tmp1(inc1+round(fl(:,1)))=fl(:,2);
tmp2(inc2+round(fd(:,1)))=fd(:,2);

fl=tmp1;
fd=tmp2;
strings=per(:);

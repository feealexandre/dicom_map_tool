function t1r_sequ_check(DATA,ROI,AddInfo,savefname)
% This Aedes plugin digs info from Siemens-T1r-prep data

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

% Grab dicom headers
p=DATA{1}.HDR.dicominfo;

% Dig TR from the first dicom header (in seconds)
TR = p(1).RepetitionTime;

% dig sWiPMemBlock:
[fl,fd,strrs]=new_parse_swip(DATA);

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


pw=WIP_PrepPulseDur;  % us, T1r PREP pulse duration
pw=pw/1000; % ms

% Number of preppulses per block
NPPB=4;  % Hard-coded in sequence: Number of Pulses Per Block

% For now, the array of pulses is as follows:
if (WIP_ArrayCheckBox==1)
    n_mlev=[0:NPPB:WIP_NumPrepPulses];
end

% create tSL-vector
tSL=pw.*n_mlev;



% generate information to display 

strs{1}=sprintf('Number of Rect Prep Pulses = %1.0f',fl(14));
strs{end+1}=sprintf('Number of Shaped Prep Pulses = %1.0f',fl(2));
strs{end+1}=sprintf('Rect Pulse Duration = %1.0f',fl(1));
strs{end+1}=sprintf('Shaped Prep Pulse Duration (T1r) = %1.0f',fl(3));
strs{end+1}=sprintf('Shaped Prep Pulse Duration (T2r) = %1.0f',fl(4));
strs{end+1}=sprintf('Shaped Prep Pulse Duration (RAFF) = %1.0f',fl(13));

strs{end+1}=sprintf('Rect Pulse Voltage Adjust = %1.0f',fd(3));
strs{end+1}=sprintf('Shaped Prep Pulse Voltage Adjust (T1r) = %1.2f',fd(1));
strs{end+1}=sprintf('Shaped Prep Pulse Voltage Adjust (T2r) = %1.2f',fd(4));
strs{end+1}=sprintf('Shaped Prep Pulse Voltage Adjust (RAFF) = %1.2f',fd(2));
strs{end+1}=sprintf('MT Offset Frequency = %1.0f',fl(8));

strs{end+1}=sprintf('Inversion pulse Applied = %1.0f    (1 = yes, 2 = no)',fl(5));
strs{end+1}=sprintf('Inversion pulse Duration = %1.0f',fl(6));
strs{end+1}=sprintf('Inversion Pulse Voltage Adjust = %1.2f',fd(5));

strs{end+1}=sprintf('Arrayed Acquisition = %1.0f     (1 = yes, 2 = no)',fl(10));
if(fl(10))
strs{end+1}=sprintf('Arrayed Size = %s',num2str(n_mlev(:)'));
end
strs{end+1}=sprintf('T2rho Applied  = %1.0f    (1 = yes, 2 = no)',fl(9));


strs{end+1}=sprintf('Spoiler Moment = %1.0f',fl(7));
strs{end+1}=sprintf('Internal TR = %1.0f',fl(11));
strs{end+1}=sprintf('Preparation Pulse Shape = %1.0f   (1=T1rHS1, 2=T1rHS4, 3=RAFF2, 4=RAFF4, 5=MT)',fl(12));


% display the information
dce_tablet(strs');





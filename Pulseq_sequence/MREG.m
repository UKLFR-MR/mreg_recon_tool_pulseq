%%This is a demo MREG sequence

%need some tools
%addpath(genpath('...\mreg_recon_tool'))

%% system limits
% maxGrad in mreg_recon_tool is 38, maxslewrate is 200, 170 not enough
dG=250e-6;
sys = mr.opts('MaxGrad', 38, 'GradUnit', 'mT/m', ...
    'MaxSlew', 170, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'adcSamplesLimit', 8192);

sys_true = mr.opts('MaxGrad', 38, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'adcSamplesLimit', 8192);

% the spoiler at the end often causes PNS problems, therefore the max slew
% rate can be set independently
sys_spoiler = mr.opts('MaxGrad', 38, 'GradUnit', 'mT/m', ...
    'MaxSlew', 60, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 10e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, 'adcSamplesLimit', 8192);

%% setup the Pulseq sequence object and the high-level parameters
seq = mr.Sequence(sys_true); %create a new sequence object
FOV = 192;          %define FOV and resolution in mm
Resolution = 3;     %mm
TR = 100e-3;         %set TR in ms

%% RF 
% due to the imperfect slice profile of the excitation pulse the slice
% thickness for the pulse is set to 0.8*slicethickness for a more
% homogeneous excitation (MZ: how is that possible? 0.8 makes it thinner, so more inhomogeneous)
flipangle = 25/180*pi;   %set flip angle in rad
sliceThickness = 150;    %set slice thickness in mm (slab excitation)
[rf, gz, gzReph] = mr.makeSincPulse(flipangle,'system',sys,'use','excitation','Duration',1.4e-3,...
    'SliceThickness',0.001*sliceThickness*0.8,'apodization',0.5,'timeBwProduct',16,'use','excitation');

%% create mreg trajectory & adc
%undersampling parameters can be varied (see stack_of_spirals.m), old
%settings are R = [3 6 2 5]. PSFs with the undersampling parameters can be
%simulated using the script stack_of_spirals_psf.m
R = [2.8687 6.0202 2.552 3.0584];

T = stack_of_spirals(R,1,1,0.001*FOV,0.001*Resolution,1,sys.maxSlew/sys.gamma, sys.maxGrad/sys.gamma*1e3); % MZ: added slew rate and gradient 
trajectStruct_export(T,'2025_Pulsec_SoS',1);% g unit T/m to mT/m, normalized
Grads_calc = -T.G'*sys.gamma;

adcSamples = length(Grads_calc)*2; % oversampling is usually 2
dW=sys.gradRasterTime/2;
% the adc should be splittable into N equal parts, each of which is aligned
% to the gradient raster. each segment however needs to have the number of
% samples divisible by 4 to be executable on siemens scanners
adcSamples=floor(adcSamples/8)*8; 
while adcSamples>0
    adcSegmentFactors=factor(round(adcSamples/2));
    assert(adcSegmentFactors(1)==2); 
    assert(adcSegmentFactors(2)==2); 
    assert(length(adcSegmentFactors)>=3); 
    adcSegments=1;
    for i=3:length(adcSegmentFactors) 
        adcSegments=adcSegments*adcSegmentFactors(i);
        adcSamplesPerSegment=adcSamples/adcSegments;
        if (adcSamplesPerSegment<=8192 && adcSegments<=128)
            break
        end
    end
    if (adcSamplesPerSegment<=8192 && adcSegments<=128)
        break
    end
    adcSamples=adcSamples-8; %
end
assert(adcSamples>0); % we could not find a suitable segmentation...
assert(adcSegments<=128);
adc = mr.makeAdc(adcSamples,'Duration',dW*adcSamples,'system',sys_true);

mregx = mr.makeArbitraryGrad('x',Grads_calc(1,:),'first',0, 'last',0,'Delay',adc.delay,'system',sys_true); % looks like the optimized MREG trajectory often exceeds the slew rate so we use an additional system with true maximum limits
mregy = mr.makeArbitraryGrad('y',Grads_calc(2,:),'first',0, 'last',0,'Delay',adc.delay,'system',sys_true);
mregz = mr.makeArbitraryGrad('z',Grads_calc(3,:),'first',0, 'last',0,'Delay',adc.delay,'system',sys_true);

%% Create fat-sat pulse 
% % not used
% % (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
B0=2.89; % 1.5 2.89 3.0
sat_ppm=-3.2;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

%% spoiler
gx_spoil=mr.makeTrapezoid('x',sys,'Area',-10*1000,'system',sys_spoiler);% spoilers cause a lot of PNS,therefore extra parameterset sys_spoiler
gy_spoil=mr.makeTrapezoid('y',sys,'Area',-10*1000,'system',sys_spoiler);
gz_spoil=mr.makeTrapezoid('z',sys,'Area',-10*1000,'system',sys_spoiler);


%trig = mr.makeDigitalOutputPulse('ext1','duration',100e-6); %'osc0', 'osc1', 'ext1' 
%% seq Block
seq = mr.Sequence(sys);
%seq.addBlock(rf_fs,gz_fs)
%seq.addBlock(gz_spoil)
seq.addBlock(rf,gz)
seq.addBlock(gzReph)

%turn off FOV shifting for readout
nopos_label = mr.makeLabel('SET', 'NOPOS', 1);
seq.addBlock(nopos_label);
seq.addBlock(mregx,mregy,mregz,adc)
%switch back on for spoiler
nopos_label_off = mr.makeLabel('SET', 'NOPOS',0);
seq.addBlock(nopos_label_off);
seq.addBlock(gx_spoil,gy_spoil,gz_spoil)
delayTR = TR - mr.calcDuration(mregx,adc)-mr.calcDuration(rf,gz)-mr.calcDuration(gzReph)-mr.calcDuration(gz_spoil);
seq.addBlock(mr.makeDelay(delayTR));

TEs = mr.calcDuration(rf, gz) + mr.calcDuration(gzReph);
TEs(:,2) = T.TE  .* 0.000001 + TEs(:,1);
seq.setDefinition('TE',TEs);% TEs is the array of TE (max. array size = 32, only first 2 matter for MREG). Unit: second
seq.plot('showblocks',1)
%% new single-function call for trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%plot trajectory if wanted
%figure; hold on;
%%plot3(ktraj(1,:),ktraj(2,:),ktraj(3,:),'b')
%plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'r'); % a 2D plot
%    %axis([-200 200 -200 200 -200 200]);
%    axis equal
%    view([0 2])



return
%% PNS 
%opportunity to simulate pns for scanner used
[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma Freiburg

%% prepare to write
[ok, error_report]=seq.checkTiming;

seq.setDefinition('FOV', 0.001*[FOV FOV sliceThickness]);
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); 
seq.setDefinition('Name', 'mreg_FOV');
seq.setDefinition('ReceiverGainHigh',1) ;
% seq.write('myMREG_FOV.seq');   % Output sequence for scanner
% %seq.install('Siemens')
% rep = seq.testReport;
% fprintf([rep{:}]);

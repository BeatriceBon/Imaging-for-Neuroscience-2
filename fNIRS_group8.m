clc
clearvars
close all

%% Loading data
%% Loading data
%setting path of useful functions  
base_path = pwd;
addpath(genpath(fullfile(base_path, 'homer2')))
addpath(genpath(fullfile(base_path, 'MNI')))
addpath(genpath(fullfile(base_path, 'iso2mesh-master')))

load('S08_walking_texting.nirs','-mat');
nCh = size(SD.MeasList,1)/2;

%% 1 - Plot 3D array configuration

figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
end
title('array configuration')
legend('sources', 'detectors', 'channels')
xlabel('x')
ylabel('y')
zlabel('z')

clear aux iCh src det

%% 2 - Compute source-detector distances and plot them with a histogram

distCh = zeros(nCh,1);
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    distCh(iCh) = sqrt(sum((src-det).^2));
end
figure;
histogram(distCh, 16) 
title('source-detectors distances histogram')
xlabel('[mm]')

clear iCh src det nbins

%% 3 - Remove noisy channels
dRange = [500 1e10];
SNRrange = 0;
remCh = removeNoisyChannels(d,dRange,SNRrange);

% Plotting array configuration highlighting bad channels in magenta
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh(iCh) == 0 % bad channels 
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'m')
    else % good channels
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
title('array configuration after having selected bad channels (magenta)')

clear iCh src det dRange SNRrange

%% 4 - Pre-process the fNIRS data

%% 4a - Convert to optical density changes
meanValue = mean(d);
dodConv = -log(abs(d)./meanValue);

clear meanValue

%% 4b - Motion correction performed with the wavelet motion correction with iqr=0.5
iqr = 0.5;
SD.MeasListAct = remCh;
dodWavelet = hmrMotionCorrectWavelet(dodConv,SD,iqr);
save ('dodWavelet.mat','dodWavelet')
% load dodWavelet.mat   % we used it to run the code faster 

clear iqr dodConv remCh

%% 4c - Band-pass filtering with cut-off frequency 0.01 and 0.5 Hz

lowerCutOff = 0.01;
higherCutOff = 0.5;
fs = 1/(t(2)-t(1)); % compute sampling frequency
dodFilt = hmrBandpassFilt(dodWavelet,fs,lowerCutOff,higherCutOff);

clear lowerCutOff higherCutOff

%% 4d - Computation of the average optical density hemodynamic response
tRange = [-2 40]; % range of timimg around stimulus to define a trial
sRange = fix(tRange*fs); % convert the time in seconds to samples
tHRF = tRange(1):1/fs:tRange(2); % time vector for the hemodynamic response (and trials)
dodAvg = zeros(length(tHRF),size(dodFilt,2),size(s,2)); % initialize the matrix that will contain our average hemodynamic response for each channel (for both wavelength) and condition
for iS = 1:size(s,2) % for each condition
    % Get the timing of stimulus presentation for that condition
    stimulusTiming = find(s(:,iS)==1); 
    % Initialize the matrix that will contain the single trial responses
    % for that condition for all channels at both wavelengths
    ytrial = zeros(length(tHRF),size(dodFilt,2),length(stimulusTiming));
    
    nTrial = 0;
    for iT = 1:length(stimulusTiming) % for each stimulus presented (for eacht trial)
        if (stimulusTiming(iT)+sRange(1))>=1 && (stimulusTiming(iT)+sRange(2))<=size(dodFilt,1) % Check that there are enough data pre and post stimulus (this is useful to check that the first stimulus is presented at least 2 seconds after the start of the acquisition and that the last stimulus has at least 18 seconds of data afterwards)
            nTrial = nTrial + 1;
            ytrial(:,:,nTrial) = dodFilt(stimulusTiming(iT)+(sRange(1):sRange(2)),:); % extract the trial from the dc data
        end
    end
    
    % Average trials (the fourth dimension of the ytrial matrix)
    dodAvg(:,:,iS) = mean(ytrial(:,:,1:nTrial),3);
    % Correct for the baseline
    for ii = 1:size(dodAvg,2) % for each channel
        foom = mean(dodAvg(1:-sRange(1),ii,iS),1); % compute baseline as average of the signal in the -2:0 seconds time range
        dodAvg(:,ii,iS) = dodAvg(:,ii,iS) - foom; % subtract the baseline from the average hemodynamic responses
    end
end

clear ii iS iT ytrial nTrial tHRF sRange foom stimulusTiming

%% 5 - Display the whole array sensitivity for the first wavelength on the volumetric GM mesh

%% Loading meshes, 10-5 positions, cranial landmarks and Jacobian matrix
load(fullfile('MNI','MNI','HeadVolumeMesh.mat'))
load(fullfile('MNI','MNI','GMSurfaceMesh.mat'))
load(fullfile('MNI','MNI','ScalpSurfaceMesh.mat'))
load(fullfile('MNI','MNI','TissueMask.mat'))

load('S08_walking_texting.jac','-mat');

clear tmp fid fileName

%% Array sensitivity on GM volume mesh for the first wavelength with every channel
HeadVolumeMesh.node(:,4) = (sum(J{1}.vol));

figure;
plotmesh(HeadVolumeMesh.node,HeadVolumeMesh.elem(HeadVolumeMesh.elem(:,5)==4,1:4))
caxis([-3 0])
colorbar
title('whole array sensitivity on the GM volume mesh, wavelength 1, all channels')

%% Array sensitivity on GM volume mesh for the first wavelength with just good channels

% Remove bad channels from Jacobian
JCropped=cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) % For each wavelength
    tmp = J{i}.vol;
    JCropped{i} = tmp(SD.MeasListAct(SD.MeasList(:,4)==i)==1,:);
end
HeadVolumeMesh_good = HeadVolumeMesh;
HeadVolumeMesh_good.node(:,4) = (sum(JCropped{1}));

figure;
plotmesh(HeadVolumeMesh_good.node,HeadVolumeMesh_good.elem(HeadVolumeMesh_good.elem(:,5)==4,1:4))
caxis([-3 0])
colorbar
title('whole array sensitivity on the GM volume mesh, wavelength 1, non-removed channels')

clear i HeadVolumeMesh_good

% the array sensitivity doesn't appear to be too different compared to the
% previous one

%% 6 - Reconstruct HbO and HbR images for both condition 1 and 2 and plot them

%% HbO and HbR images reconstruction
% Compute inverse of Jacobian
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our
% data
datarecon = -dodAvg;

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(SD.Lambda); % n of wavelengths
nCond = size(s,2); % number of condition

% Initialize final results matrices for condition 1 (texting on right hand)
hbo.vol = zeros(nFrames,nNodeVol,nCond);
hbr.vol = zeros(nFrames,nNodeVol,nCond);
hbo.gm = zeros(nFrames,nNodeGM,nCond);
hbr.gm = zeros(nFrames,nNodeGM,nCond);


% Obtain specific absorption coefficients
Eall = zeros(nWavs,2);
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall(i,:) = Etmp./1e7; %This will be nWavs x 2;
end

% Perform reconstruction for each frame
    
for cond = 1:nCond
    
    % For each frame
    for frame = 1:nFrames
        
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & SD.MeasListAct==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
        
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);        
        
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
        
        % Book-keeping and saving
        hbo.vol(frame,:,cond) = hbo_tmpVol;
        hbr.vol(frame,:,cond) = hbr_tmpVol;
        hbo.gm(frame,:,cond) = hbo_tmpGM;
        hbr.gm(frame,:,cond) = hbr_tmpGM;
        
    end
end

clear lambda1 invJ i Jtmp JJT S nNodeGM wavelengths Etmp datarecon datarecon_regr Eall nNodeVol
clear cond frame muaImageAll wav dataTmp invJtmp tmp hbo_tmpVol hbr_tmpVol hbo_tmpGM hbr_tmpGM vol2gm

%% Plotting HbO and HbR reconstructed images for both conditions at time points 0s, 10s and 18s
tRecon = [0 10 18];
baseline = abs(tRange(1)); % two seconds of baseline
sRecon = fix(tRecon*fs)+fix(baseline*fs); % Convert to samples
load greyJet

for iT = 1:length(sRecon) % for each time point
    for iCond = 1:nCond % for each condition
        
        figure
        subplot(121)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbo.gm(sRecon(iT),:,iCond);
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        caxis([-0.2 0.2]) % Set the limit of the colorbar
        view([0 90]) % Set the view angle
        title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet) % set the loaded colormap
        hb = colorbar;
        hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
        axis off % remove axis

        subplot(122)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbr.gm(sRecon(iT),:,iCond);
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        view([0 90])
        caxis([-0.2 0.2])
        title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet)
        hb = colorbar;
        hb.Label.String = {'\DeltaHbR [\muM]'};
        axis off

    end
end

hbo_vol=hbo.vol(sRecon,:,:);
hbr_vol=hbr.vol(sRecon,:,:);
hbo_gm=hbo.vol(sRecon,:,:);
hbr_gm=hbr.vol(sRecon,:,:);

clear baseline iCond iT hb greyJet

%% 7 - Reconstruct HbO and HbR images again but removing channels with SD distance >=30mm

%% Removing channels with SD distance >=30mm and displaying resulting array sensitivity (not needed)
longDist=distCh<30;
longDist(nCh+1:2*nCh)=longDist(1:nCh);
remCh_new=SD.MeasListAct.*double(longDist);

JCropped=cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) % For each wavelength
    tmp = J{i}.vol;
    JCropped{i} = tmp(remCh_new(SD.MeasList(:,4)==i)==1,:);
end
HeadVolumeMesh_good = HeadVolumeMesh;
HeadVolumeMesh_good.node(:,4) = (sum(JCropped{1}));

% Plotting new array configuration, in magenta the removed channels
figure;
plot3(SD.SrcPos(:,1),SD.SrcPos(:,2),SD.SrcPos(:,3),'.r','MarkerSize',10)
hold on;
plot3(SD.DetPos(:,1),SD.DetPos(:,2),SD.DetPos(:,3),'.b','MarkerSize',10)
for iCh = 1:nCh
    src = SD.SrcPos(SD.MeasList(iCh,1),:);
    det = SD.DetPos(SD.MeasList(iCh,2),:);
    if remCh_new(iCh) == 0 % bad channels 
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'m')
    else % good channels
        plot3([src(1) det(1)],[src(2) det(2)],[src(3) det(3)],'g')
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
title('array configuration after having selected bad channels (magenta)')

% New array sensitivity, seems to be much worse compared to the one we had
% before, which makes sense because the spatial resolution is smaller
% than what we needed
figure;
plotmesh(HeadVolumeMesh_good.node,HeadVolumeMesh_good.elem(HeadVolumeMesh_good.elem(:,5)==4,1:4))
caxis([-3 0])
colorbar
title('whole array sensitivity on the GM volume mesh, wavelength 1, non-removed channels')

clear i HeadVolumeMesh_good src det distCh longDist


%% HbO and HbR image reconstruction
% Compute inverse of Jacobian
lambda1 = 0.1;
invJ = cell(length(SD.Lambda),1);
for i = 1:length(SD.Lambda) %for each Jacobian
    Jtmp = JCropped{i};
    JJT = Jtmp*Jtmp';
    S=svd(JJT);
    invJ{i} = Jtmp'/(JJT + eye(length(JJT))*(lambda1*max(S)));
end

% Data to reconstruct are optical density changes compared to a baseline.
% In our case the baseline is 0, therefore we want to reconstruct 0-our
% data.
datarecon = -dodAvg; 

% Inizialize matrices and load useful stuff
nNodeVol = size(HeadVolumeMesh.node,1);  %The node count of the volume mesh
nNodeGM = size(GMSurfaceMesh.node,1); %The node count of the GM mesh
nFrames = size(datarecon,1); % Number of samples to reconstruct
load('vol2gm.mat')
wavelengths = SD.Lambda; % wavelengths of the system
nWavs = length(SD.Lambda); % n of wavelengths
nCond = size(s,2); % number of condition

% Initialize final results matrices for condition 1 (texting on right hand)
hbo_close.vol = zeros(nFrames,nNodeVol,nCond);
hbr_close.vol = zeros(nFrames,nNodeVol,nCond);
hbo_close.gm = zeros(nFrames,nNodeGM,nCond);
hbr_close.gm = zeros(nFrames,nNodeGM,nCond);


% Obtain specific absorption coefficients
Eall = zeros(nWavs,2);
for i = 1:nWavs
    Etmp = GetExtinctions(wavelengths(i));
    Etmp = Etmp(1:2); %HbO and HbR only
    Eall(i,:) = Etmp./1e7; %This will be nWavs x 2;
end

% Perform reconstruction for each frame
    
for cond = 1:nCond
    
    % For each frame
    for frame = 1:nFrames
        
        % Reconstruct absorption changes
        muaImageAll = zeros(nWavs,nNodeVol);
        for wav = 1:nWavs
            dataTmp = squeeze(datarecon(frame,SD.MeasList(:,4)==wav & remCh_new==1,cond));
            invJtmp = invJ{wav};
            tmp = invJtmp * dataTmp';
            muaImageAll(wav,:) = tmp; %This will be nWavs * nNode
        end
        
        % Convert to concentration changes
        hbo_tmpVol = (Eall(2,2)*muaImageAll(1,:) - Eall(1,2)*muaImageAll(2,:))/(Eall(1,1)*Eall(2,2)-Eall(1,2)*Eall(2,1));
        hbr_tmpVol = (muaImageAll(2,:)-Eall(1,2)*hbo_tmpVol)/Eall(2,2);        
        
        % Map to GM surface mesh
        hbo_tmpGM = (vol2gm*hbo_tmpVol');
        hbr_tmpGM = (vol2gm*hbr_tmpVol');
        
        % Book-keeping and saving
        hbo_close.vol(frame,:,cond) = hbo_tmpVol;
        hbr_close.vol(frame,:,cond) = hbr_tmpVol;
        hbo_close.gm(frame,:,cond) = hbo_tmpGM;
        hbr_close.gm(frame,:,cond) = hbr_tmpGM;
        
    end
end

clear lambda1 invJ i Jtmp JJT S nNodeGM wavelengths Etmp datarecon datarecon_regr Eall nNodeVol
clear cond frame muaImageAll wav dataTmp invJtmp tmp hbo_tmpVol hbr_tmpVol hbo_tmpGM hbr_tmpGM vol2gm

%% Plotting HbO and HbR reconstructed images for both conditions at time points 0s, 10s and 18s with new array
tRecon = [0 10 18];
baseline = abs(tRange(1)); % two seconds of baseline
sRecon = fix(tRecon*fs)+fix(baseline*fs); % Convert to samples
load greyJet

for iT = 1:length(sRecon) % for each time point
    for iCond = 1:nCond % for each condition
        
        figure
        subplot(121)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbo_close.gm(sRecon(iT),:,iCond);
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        caxis([-0.2 0.2]) % Set the limit of the colorbar
        view([0 90]) % Set the view angle
        title(['HbO cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet) % set the loaded colormap
        hb = colorbar;
        hb.Label.String = {'\DeltaHbO [\muM]'}; % assign label to colorbar
        axis off % remove axis

        subplot(122)
        % Assign image to fourth column of node
        GMSurfaceMesh.node(:,4) = hbr_close.gm(sRecon(iT),:,iCond);
        plotmesh(GMSurfaceMesh.node,GMSurfaceMesh.face)
        view([0 90])
        caxis([-0.2 0.2])
        title(['HbR cond ' num2str(iCond) ', t = ' num2str(tRecon(iT)) ' s'])
        colormap(greyJet)
        hb = colorbar;
        hb.Label.String = {'\DeltaHbR [\muM]'};
        axis off

    end
end

hbo_close_vol=hbo_close.vol(sRecon,:,:);
hbr_close_vol=hbr_close.vol(sRecon,:,:);
hbo_close_gm=hbo_close.vol(sRecon,:,:);
hbr_close_gm=hbr_close.vol(sRecon,:,:);
remCh=SD.MeasListAct;

clear baseline iCond iT hb greyJet


%% 8 - Compare the reconstructions from points 6 and 7
% The reconstruction made without removing the channels with distance >=30mm
% is the more accurate one since, as we can also see from the array
% sensitivity, it has a much better resolution of the brain activity.
% Furthermore, if HbO, for example, increases HbR is supposed to decrease and
% viceversa, but this doesn't happen in the case of the array in which we removed
% the long distance channels and this shows that the images given by this array
% are not as good as the one given by the one in point 6.

save Results/results_fNIRS.mat remCh hbo_vol hbr_vol hbo_gm hbr_gm hbo_close_vol hbr_close_vol hbo_close_gm hbr_close_gm



function [divergencevectorinfrequency, divergencevectorintime ] = KLDistanceforieeg()
%KLDistanceforieeg calculates the Kullback Leibier and other distances.
% This function calls to KLDistanceforieeg_onepatient(patientid)

%prepare loop to call single patients
%eegpatientl = { 'TWH027','TWH024','TWH028','TWH030', 'TWH031','TWH033','TWH034', 'TWH037','TWH038','TWH042','TWH043'};
global globalFsDir;
[globalFsDir] = loadglobalFsDir();
eegpatientl = {'TWH030','TWH031','TWH033','TWH034','TWH037','TWH038','TWH042','TWH043'};
eegconditions = {'EC_PRE' 'EC_POST'; 'EC_PRE' 'HYP'};% 'EC_PRE' 'EO_POST'};
for indpat=1:length(eegpatientl)
    patientid = eegpatientl(indpat);
    for indcond=1:length(eegconditions)
        conditionid = eegconditions(indcond,:);
        prelabel = conditionid(1);prelabel = prelabel{1};
        postlabel = conditionid(2);postlabel = postlabel{1};
        labels = {prelabel, postlabel};
        %calculate statistical distances(KL, KS, EMD)
        fprintf('Calculating Statistical Distances for patient in Frequency domain %s between %s and %s \n',patientid{1},prelabel, postlabel)
        %divergencevectorinfrequency(indpat,indcond) = KLDistanceinfrequencydom_onepatient(patientid{1},prelabel, postlabel);
        fprintf('Calculating Statistical Distances for patient in Time domain %s between %s and %s \n',patientid{1},prelabel, postlabel);
        divergencevectorintime(indpat,indcond) = KLDistanceintimedom_onepatient(patientid{1},prelabel, postlabel);
    end
end
%plotKLdistancesinfreq(divergencevectorinfrequency, labels);
plotKLdistancesintime(divergencevectorintime, eegconditions);
end
function divergence = KLDistanceintimedom_onepatient(patientid, prelabel, postlabel)
%KLDistanceintimedom_onepatient calculates KL, KS and EMD distances between
%prelabel and postlabel  in time domain from the EEG power time series
[myfullname, EEG_pre, channel_labels, eegdate, eegsession] = initialize_EEG_variables(patientid,prelabel);
[myfullname, EEG_post, channel_labels, eegdate, eegsession] = initialize_EEG_variables(patientid,postlabel);
nb_channels = length(channel_labels) -1;
divergence.patientid = patientid;
minnbcols = min(size(EEG_post.data,2),size(EEG_pre.data,2));
divergence.channel_labels = channel_labels(2:end);
psize = length(EEG_pre.data);
qsize = length(EEG_post.data);
psize = min(psize,qsize);
divergence.kl_power_ts = mean(KLDiv(abs(EEG_pre.data(:,(1:psize))),abs(EEG_post.data(:,(1:psize)))));
% kolmogorov Smirnoff
ks_total_distance = 0;
ks_total_p = 0;
for i=1:nb_channels
    [divergence.ks_h_ts,divergence.ks_p_ts,divergence.ks_distance_ts] = kstest2(EEG_pre.data(i,(1:psize)),EEG_post.data(i,(1:psize)));
    ks_total_distance =ks_total_distance + divergence.ks_distance_ts;
    ks_total_p = ks_total_p + divergence.ks_p_ts;
end
%calculate the normalized p and distance
divergence.ks_p_ts = ks_total_p/nb_channels;
divergence.ks_distance_ts = ks_total_distance/nb_channels;
%EMD distance
 weight_emd=repmat(1/nb_channels,1,nb_channels);
 fprintf('.....Take a seat, calculating the EMD of power time series in patient:%s \n', patientid)
 [x, divergence.emd_power_ts] =emd(double(EEG_pre.data(2:end,1:minnbcols)),double(EEG_post.data(2:end,1:minnbcols)),weight_emd,weight_emd,@gdf);
 fprintf('The EMD for the entire power time series between %s and %s in patient:%s is : %.2f \n',prelabel,postlabel, patientid,divergence.emd_power_ts)
end
function [] = plotKLdistancesintime(vectorpatientdistances, labels)
%plotKLdistances  plot distances in time domain

nb_pats = size(vectorpatientdistances,1);
nb_pairofconditions = size(vectorpatientdistances,2);
legend_text = cell(1,nb_pairofconditions);
legend_text{1}= 'EC PRE-EC POST';legend_text{2}= 'ECPRE-HYP';legend_text{3}= 'EC PRE-EO PRE';
%x=[1:1:nb_pats]';
%plot the average per channel for all patients
maxpower_kl = 0; maxpower_ks = maxpower_kl;
ykl_all = []; yks_all= [];yemd_all = [];
for i=1:nb_pats
    ykl = [];yks = []; yemd= [];
    for j=1:nb_pairofconditions
        patdist_kl = vectorpatientdistances(i,j).kl_power_ts;
        patdist_ks =  vectorpatientdistances(i,j).ks_distance_ts;
        patdist_emd = vectorpatientdistances(i,j).emd_power_ts;
        ykl = [ykl; patdist_kl];
        yks = [yks; patdist_ks];
        yemd = [yemd; patdist_emd];
    end
    ykl_all = [ykl_all; ykl];yks_all = [yks_all; yks];yemd_all = [yemd_all; yemd];
    %plot single patient all pairs of conditions
    hdist_kl = figure;
    bar(ykl);
    grid on
    set(gca, 'YLim', [0 1.5*max(ykl)],'xticklabel', legend_text);
    
    %xticklabel_rotate([],45,[],'Fontsize',6);
    xlabel('Conditions Pair');ylabel('KL Distance')
    msgtitle = sprintf('Patient = %s, KL Distance of power time series, mean for all channels', vectorpatientdistances(i,j).patientid);
    title(msgtitle);
    
    hdist2 = figure;
    bar(yks);
    grid on;
    set(gca, 'YLim', [0 1.5*max(yks)],'xticklabel', legend_text);
    %xticklabel_rotate([],45,[],'Fontsize',6);
    xlabel('Conditions Pair');ylabel('KS Distance')
    msgtitle = sprintf('Patient = %s, KS Distance of power time series, mean for all channels', vectorpatientdistances(i,j).patientid);
    title(msgtitle);
    %EMD one patient
    hdist3 = figure;
    bar(yemd);
    grid on;
    set(gca, 'YLim', [0 1.5*max(yemd)],'xticklabel', legend_text);
    %xticklabel_rotate([],45,[],'Fontsize',6);
    xlabel('Conditions Pair');ylabel('EMD Distance')
    msgtitle = sprintf('Patient = %s, EMD Distance of power time series, mean for all channels', vectorpatientdistances(i,j).patientid);
    title(msgtitle);
end

hdistkla = figure;
bar(reshape(ykl_all,[nb_pairofconditions, nb_pats])');
set(gca, 'YLim', [0 1.5*max(ykl_all)],'XTickLabel',{vectorpatientdistances(1:end).patientid});
%xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients');ylabel('KL Distance')
msgtitle = sprintf('KL Distance all patients');
title(msgtitle);
legend(legend_text{1},legend_text{2},legend_text{3});
grid on;

%Chart for KS distance
hdistksa = figure;
bar(reshape(yks_all,[nb_pairofconditions, nb_pats])');
set(gca, 'YLim', [0 1.5*max(yks_all)],'XTickLabel',{vectorpatientdistances(1:end).patientid});
%xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients');ylabel('KS Distance')
msgtitle = sprintf('KS Distance  all patients');
title(msgtitle);
legend(legend_text{1},legend_text{2},legend_text{3});
grid on;
hdistemda = figure;
bar(reshape(yemd_all,[nb_pairofconditions, nb_pats])');
set(gca, 'YLim', [0 1.5*max(yemd_all)],'XTickLabel',{vectorpatientdistances(1:end).patientid});
%xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients'); ylabel('EMD Distance');
msgtitle = sprintf('EMD Distance  all patients');
title(msgtitle);
legend(legend_text{1},legend_text{2},legend_text{3});
grid on;
end



function divergence = KLDistanceinfrequencydom_onepatient(patientid, prelabel, postlabel)
% KLDistanceforieeg_onepatient calculates KL, KS and EMD distances between
%prelabel and postlabel for several variables in the mat files fft_EC_POST_TWH043_05042016_s1
% Example:
% divergence = KLDistanceforieeg_onepatient('TWH043','EC','EC')
% requested_frequences_power_bnds, requested_frequences_power_bnds(),frqperbad
% requested_frequences_power_bnds(channels, bands)
% requested_frequences_power_bnds(2,:) = 5.0972    0.7405    0.0284    0.0205    0.0057
% channel 2 has 5.0972 of power in absolute value. To know how much in each
% band each channel picks up:percentlistoffrqperband
% percentlistoffrqperband in delta for channel_2 = 5.0972/sum(5.0972 ,0.7405,0.0284,0.0205,0.0057)
% each band
% global globalFsDir;
% [globalFsDir] = loadglobalFsDir();
patdirpatient = fullfile(globalFsDir,patientid, 'data\figures');
%[pre,post] = getpreandpostmatfiles(patientid);
[pre, post] = gettwosampleseries(patientid, prelabel, postlabel);
%Build distribution q
fh_q = fullfile(patdirpatient,post);
q_vofi  = load(fh_q);

%Build distribution p
fh_p = fullfile(patdirpatient,pre);
p_vofi  = load(fh_p);
% Build divergence vector
%shortest time epoch necessary to calculate KL because if for equal length
%vectors
nb_channels = length(p_vofi.requested_frequences_power_bnds);
nb_bands = length(p_vofi.frqperband);
mintimets = min(length(p_vofi.power_fft), length(q_vofi.power_fft));
divergence.patientid = patientid;
divergence.channel_labels = p_vofi.channel_labels;
divergence.kl_power_ts = KLDiv(p_vofi.power_fft(1:mintimets),q_vofi.power_fft(1:mintimets));
divergence.kl_div_frqperband = KLDiv(p_vofi.frqperband',q_vofi.frqperband');
divergence.kl_requested_frequences_power_bnds_pc = KLDiv(p_vofi.requested_frequences_power_bnds, q_vofi.requested_frequences_power_bnds);

for i=1:nb_channels
    for j=1:nb_bands
        p_availablepowerinaband_pc(i,j) = p_vofi.requested_frequences_power_bnds(i,j)/sum(p_vofi.requested_frequences_power_bnds(:,j));
        q_availablepowerinaband_pc(i,j) = q_vofi.requested_frequences_power_bnds(i,j)/sum(q_vofi.requested_frequences_power_bnds(:,j));
    end
end
divergence.kl_availablepowerinaband_pc = KLDiv(p_availablepowerinaband_pc,q_availablepowerinaband_pc);
% kolmogorov Smirnoff
ks_total_distance = 0;
ks_total_p = 0;
for i=1:nb_channels
    [divergence.ks_h_ts,divergence.ks_p_ts,divergence.ks_distance_ts] = kstest2(p_vofi.power_fft(i,:),q_vofi.power_fft(i,:));
    ks_total_distance =ks_total_distance + divergence.ks_distance_ts;
    ks_total_p = ks_total_p + divergence.ks_p_ts;
end
%calculate the normalized p and distance
divergence.ks_p_ts = ks_total_p/nb_channels;
divergence.ks_distance_ts = ks_total_distance/nb_channels;
%convert matrix nb_channelsxfreq bands into the 1 (mean)XFreq
p_power_perbandmean_v = [mean(p_vofi.requested_frequences_power_bnds(:,1)),mean(p_vofi.requested_frequences_power_bnds(:,2)),mean(p_vofi.requested_frequences_power_bnds(:,3)),mean(p_vofi.requested_frequences_power_bnds(:,4)),mean(p_vofi.requested_frequences_power_bnds(:,5))];
q_power_perbandmean_v = [mean(q_vofi.requested_frequences_power_bnds(:,1)),mean(q_vofi.requested_frequences_power_bnds(:,2)),mean(q_vofi.requested_frequences_power_bnds(:,3)),mean(q_vofi.requested_frequences_power_bnds(:,4)),mean(q_vofi.requested_frequences_power_bnds(:,5))];
[divergence.ks_h_meanpower_bnds_pc,divergence.ks_p_meanfrequences_power_bnds_pc,divergence.ks_dist_meanfrequences_power_bnds_pc] = kstest2(p_power_perbandmean_v, q_power_perbandmean_v);
% the emd distance.
weight_emd=repmat(1/nb_channels,1,nb_channels);
% x is the flow that shows the flow that minimizes the flow
%cost for the transportation problem converges and fval is the
% value of this flow.
fprintf('.....Take a seat, calculating the EMD of power in frequency domain in patient:%s \n', patientid)
[x, divergence.emd_power_ts] = emd(p_vofi.power_fft,q_vofi.power_fft,weight_emd,weight_emd,@gdf);
fprintf('The EMD for the entire power time series between %s and %s in patient:%s is : %.2f \n',prelabel,postlabel, patientid,divergence.emd_power_ts)
end

function [] = plotKLdistancesinfreq(vectorpatientdistances, labels)
% plotKLdistances plot KL distance according to the vectorpatientdistances
%1 KL, 2 KS 3 EMD
%1. Charts for KL distance
hchan = figure;
% hpat = figure;
nb_pats = length(vectorpatientdistances);
nb_bands = length(vectorpatientdistances(1).kl_requested_frequences_power_bnds_pc);
subsetchannels = 0;
if subsetchannels > 0
    hchanregio = figure;
end
%plot the average per channel for all patients
maxpower = 0; nb_channels = [];
for i=1:length(vectorpatientdistances)
    nb_channels(i) = length(vectorpatientdistances(i).channel_labels) -1;
    if maxpower < sum(vectorpatientdistances(i).kl_availablepowerinaband_pc)/nb_channels(i)
        maxpower =  sum(vectorpatientdistances(i).kl_availablepowerinaband_pc)/nb_channels(i);
    end
end
ypc = []; ypb = []; localel =[];maxpowerdivband =0; ypc_ks = [];ypc_emd = [];
for i=1:length(vectorpatientdistances)
    ypc = [ypc; sum(vectorpatientdistances(i).kl_availablepowerinaband_pc)/nb_channels(i)];
    ypc_ks = [ypc_ks;vectorpatientdistances(i).ks_distance_ts ];
    ypc_emd = [ypc_emd;vectorpatientdistances(i).emd_power_ts ];
    %     localel = vectorpatientdistances(i).kl_requested_frequences_power_bnds_pc;
    %     if maxpowerdivband < max(localel)
    %            maxpowerdivband = max(localel);
    %     end
    %     ypb = [ypb;localel'];
    %     localel =[];
end
bar(ypc);
set(gca, 'YLim', [0 1.5*maxpower],'XTickLabel',{vectorpatientdistances(1:end).patientid});
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients');ylabel('KL Distance')
msgtitle = sprintf('KL Distance %s vs %s of power distribution in frequency domain', labels{1}, labels{2});
title(msgtitle);

% plot frequency band distance
% hband = figure;
% bar(ypb');
% set(gca, 'YLim', [0 1.5*maxpowerdivband],'XTick',[1:5],'XTickLabel',{'delta' 'theta' 'alpha' 'beta' 'gamma'});
% xlabel('Bands'); ylabel('KL Distance per band')
% msgtitle = sprintf('KL Distance Eyes Closed PRE-POST of power distribution per band');
% legend({vectorpatientdistances(1:end).patientid})
% title(msgtitle);

%2. KS distance
hks = figure;
bar(ypc_ks);
set(gca, 'YLim', [0 1.5*max(ypc_ks)],'XTickLabel',{vectorpatientdistances(1:end).patientid});
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients');ylabel('KL Distance')
msgtitle = sprintf('KS Distance %s vs %s of power distribution in frequency domain',labels{1}, labels{2});
title(msgtitle);

%3. EMD distance
hemd = figure;
bar(ypc_emd);
set(gca, 'YLim', [0 1.5*max(ypc_emd)],'XTickLabel',{vectorpatientdistances(1:end).patientid});
xticklabel_rotate([],45,[],'Fontsize',6);
xlabel('Patients');ylabel('EMD Distance')
msgtitle = sprintf('EMD Distance %s vs %s of power distribution in frequency domain',labels{1}, labels{2});
title(msgtitle);
% print out divergence structure
%disp(divergence);
end
function [pre_series, post_series] = gettwosampleseries(patientid, prelabel, postlabel)
%gettwosampleseries
% global globalFsDir;
% [globalFsDir] = loadglobalFsDir();
% patdirpatient = fullfile(globalFsDir,patientid, 'data\figures');
[pre_series,post_series] = getpreandpostmatfiles(patientid, prelabel, postlabel);
end
function [pre,post] = getpreandpostmatfiles(patientid, pre_label, post_label)
% getpreandpostmatfiles return the mat file that contains the fft results
% for power
% Example:
% [pre,post] = getpreandpostmatfiles(TWH043)
%  pre = 'fft_EC_PRE_TWH043_05042016_s1.mat'; post = 'fft_EC_POST_TWH043_05042016_s1.mat';
if strcmp(patientid, 'TWH043') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH043_05042016_s1.mat';
        post = 'fft_EC_POST_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH043_05042016_s1.mat';
        post = 'fft_HYP_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH043_05042016_s1.mat';
        post = 'fft_HYP_TWH043_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH043_05042016_s1.mat';
        post = 'fft_EC_POST_TWH043_05042016_s1.mat';
    end
elseif strcmp(patientid, 'TWH042') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH042_05042016_s1.mat';
        post = 'fft_EC_POST_TWH042_05042016_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH042_05042016_s1.mat';
        post = 'fft_HYP_TWH042_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH042_05042016_s1.mat';
        post = 'fft_HYP_TWH042_05042016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH042_05042016_s1.mat';
        post = 'fft_EC_POST_TWH042_05042016_s1.mat';
    end
elseif strcmp(patientid, 'TWH038') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH038_03082016_s1.mat';
        post = 'fft_EC_POST_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH038_03082016_s1.mat';
        post = 'fft_HYP_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH038_03082016_s1.mat';
        post = 'fft_HYP_TWH038_03082016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH038_03082016_s1.mat';
        post = 'fft_EC_POST_TWH038_03082016_s1.mat';
    end
elseif strcmp(patientid, 'TWH037') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH037_03142016_s1.mat';
        post = 'fft_EC_POST_TWH037_03142016_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH037_03142016_s1.mat';
        post = 'fft_HYP_TWH037_03142016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH037_03142016_s1.mat';
        post = 'fft_HYP_TWH037_03142016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH037_03142016_s1.mat';
        post = 'fft_EC_POST_TWH037_03142016_s1.mat';
    end
elseif strcmp(patientid, 'TWH034') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH034_02092016_s2.mat';
        post = 'fft_EC_POST_TWH034_02092016_s2.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH034_02092016_s2.mat';
        post = 'fft_HYP_TWH034_02092016_s2.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH034_02092016_s2.mat';
        post = 'fft_HYP_TWH034_02092016_s2.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH034_02092016_s2.mat';
        post = 'fft_EC_POST_TWH034_02092016_s2.mat';
    end
elseif strcmp(patientid, 'TWH033') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH033_02032016_s1.mat';
        post = 'fft_EC_POST_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH033_02032016_s1.mat';
        post = 'fft_HYP_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH033_02032016_s1.mat';
        post = 'fft_HYP_TWH033_02032016_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH033_02032016_s1.mat';
        post = 'fft_EC_POST_TWH033_02032016_s1.mat';
    end
elseif strcmp(patientid, 'TWH031') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH031_12012015_s1.mat';
        post = 'fft_EC_POST_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH031_12012015_s1.mat';
        post = 'fft_HYP_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH031_12012015_s1.mat';
        post = 'fft_HYP_TWH031_12012015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH031_12012015_s1.mat';
        post = 'fft_EC_POST_TWH031_12012015_s1.mat';
    end
elseif strcmp(patientid, 'TWH030') ==1
    if strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EC_PRE_TWH030_11172015_s1.mat';
        post = 'fft_EC_POST_TWH030_11172015_s1.mat';
    elseif strcmp(pre_label, 'EC_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EC_PRE_TWH030_11172015_s1.mat';
        post = 'fft_HYP_TWH030_11172015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'HYP') > 1
        pre = 'fft_EO_PRE_TWH030_11172015_s1.mat';
        post = 'fft_HYP_TWH030_11172015_s1.mat';
    elseif strcmp(pre_label, 'EO_PRE') + strcmp(post_label, 'EC_POST') > 1
        pre = 'fft_EO_PRE_TWH030_11172015_s1.mat';
        post = 'fft_EC_POST_TWH030_11172015_s1.mat';
    end
end
end
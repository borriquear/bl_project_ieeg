function [myfullname, EEG, channel_labels] =  initialize_EEG_variables()
clear all
% initialize_EEG_variables  load mat file and initialize EEG lab variables .
%   [myfullname, EEG_study, channel_labels] =
%   initialize_EEG_variables() load matfile, save the EEG object as
%   EEG_study. Note this is hard-code: EEG_study = EEG_cut_BL_HYP
%   channel_labels labels of the channels
%   myfullname contains information of the condition,session, patient and date 
%
%Load matfile 

mydir = 'D:\BIAL PROJECT\patients\TWH027\data'
mydir = 'D:\BIAL PROJECT\patients\TWH028\data'

mydir = 'D:\BIAL PROJECT\patients\TWH033\data'
mydir = 'D:\BIAL PROJECT\patients\TWH034\data'
mydir = 'D:\BIAL PROJECT\patients\TWH024\data'
mydir = 'D:\BIAL PROJECT\patients\TWH030\data'
%myfile = 'EEG_cut_BL_EO_POST_mj34_02092016_s2.mat'

myfile = 'EEG_cut_BL_HYP_bs27_10222015_s2.mat'
myfile = 'EEG_cut_BL_HYP_cj28_10202015_s1'

myfile = 'EEG_cut_BL_HYP_nk33_02032016_s1.mat'
myfile = 'EEG_cut_BL_HYP_mj34_02092016_s2'
myfile = 'EEG_cut_BL_HYP_fo24_09192015_s1.mat'
myfile = 'EEG_cut_BL_HYP_was30_11172015_s1.mat'

myfullname = fullfile(mydir, myfile);
disp('Loading mat file...')
EEG = load(myfullname);
if isstruct(EEG) == 1
    if (strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH033\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH034\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH024\data')) == 1
        EEG=EEG.EEG
    elseif (strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH030\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH027\data') + strcmp(mydir, 'D:\BIAL PROJECT\patients\TWH028\data')) == 1
        EEG = EEG.EEG_cut_BL_HYP
        %EEG = EEG.EEG_cut_BL_EC_PRE
        %EEG = EEG.EEG_cut_BL_EO_PRE
        %EEG = EEG.EEG_cut_BL_EC_POST
    end
end
disp(['File ' myfile ' loaded!' ])
channel_labels = {EEG.chanlocs.labels};
disp([ 'Displaying the label of all the channels....' ])
initchan = 2 % channel 1 is the null Event channel
for i=initchan:EEG.nbchan
    disp(['Index ' i ' is the channel' channel_labels(i)])
end

end


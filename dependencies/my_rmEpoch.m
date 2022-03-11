function [output,rm_idx] = my_rmEpoch(EEG)
EEG = pop_jointprob(EEG,1,1,5,5,0,0,0,[],0);
EEG = pop_eegthresh(EEG,1,1,-100,100,-0.5,0.998,2,0);
EEG = pop_rejtrend(EEG,1,1,750,50,0.3,2,0);
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
output = pop_rejepoch(EEG,EEG.reject.rejglobal,0);
rm_idx = EEG.reject.rejglobal;
end
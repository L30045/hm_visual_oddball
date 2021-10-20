function EEG = insert_data(EEG, data, ch_name)
EEG.nbchan = EEG.nbchan +1;
EEG.data = [EEG.data;data];
EEG.chanlocs(end+1).labels = ch_name;
end
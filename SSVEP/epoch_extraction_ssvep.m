%% Script for epoching SSVEP data
filepath = 'C:\Users\Yuan\OneDrive\Desktop\sub1163\';
filename = '2004_SSVEP_1163_condition.xdf';

rmbase_flag = true;
rmtrial_flag = false;

epoch_struct = epoch_ssvep([filepath,filename],rmbase_flag,rmtrial_flag);

%% FFT
test_epoch = epoch_struct.ring1.stim.up;
tarCh = {'O1','O2','Oz'};
ev_name = test_epoch.event(find(cellfun(@(x) ~isempty(regexp(x,'Frequency','once')), {test_epoch.event.type}),1)).type;
ev_name = regexp(ev_name,'\;','split');
ev_name = ev_name{4};
tarFreq = str2double(ev_name(regexp(ev_name,'\:')+2:end));
plt_ch_idx = ismember({test_epoch.chanlocs.labels},tarCh);
tar_data = test_epoch.data(ismember({test_epoch.chanlocs.labels},tarCh),:,:);
% FFT
Fs = test_epoch.srate;
plt_t = test_epoch.times;
plt_f = 1:20; % upper Hz
plt_fft = zeros(length(plt_f),size(tar_data,3)); % freq by trials

for t_i = 1:size(tar_data,3)
    tmp_fft = fft(tar_data(:,:,t_i)');
    L = size(tmp_fft,1);
    f = Fs*(0:(L/2))/L;
    P2 = abs(tmp_fft/L);
    P1 = P2(1:L/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);
    f_idx = ismember(f,plt_f);
    plt_fft(:,t_i) = mean(P1(f_idx,:),2);
end

figure
plot(plt_f, mean(plt_fft,2),'linewidth',3,'displayname',sprintf('Target Freq: %d',tarFreq))
legend
grid on
xlabel('Frequency (Hz)')
ylabel('Amplitude (V^2)')
set(gca,'fontsize',10)


%% Behavior
ang_lib = cell(4, 2); % direction * ring
angDiff_lib = cell(4, 2); % direction * ring
dist_lib = cell(4, 2); % direction * ring

for cond_i = 1:2    
    upLoc = epoch_struct.event_time.upLoc;
    downLoc = epoch_struct.event_time.downLoc;
    leftLoc = epoch_struct.event_time.leftLoc;
    rightLoc = epoch_struct.event_time.rightLoc;
    tar_lib = [upLoc;downLoc;leftLoc;rightLoc];
%         ev_list = fieldnames(epoch_struct);
    ev_list = {'std_epoch','dev_epoch','gip_std','gip_dev'};
    ev_direct = {[epoch_struct.event_time.std_up;epoch_struct.event_time.std_down;...
                  epoch_struct.event_time.std_left; epoch_struct.event_time.std_right];...
                 [epoch_struct.event_time.dev_up;epoch_struct.event_time.dev_down;...
                  epoch_struct.event_time.dev_left; epoch_struct.event_time.dev_right]};
    dir_idx = [1,2,1,2];

    for e_i = 1:4
        tar_epoch = epoch_struct.(ev_list{e_i});
        nbchan = find(ismember({tar_epoch.chanlocs.labels},'HeadLoc_x'));
        tar_direct = ev_direct{dir_idx(e_i)};
        if e_i == 3
            s_t = isnan(epoch_struct.event_time.diff_gip_std);
        elseif e_i == 4
            s_t = isnan(epoch_struct.event_time.diff_gip_dev);
        else
            s_t = false(1,size(tar_direct,2));
        end
        tar_direct(:,s_t) = [];
        headLoc = tar_epoch.data(nbchan:nbchan+2,:,:);
        headDirect = tar_epoch.data(nbchan+3:nbchan+5,:,:);
        GIP = tar_epoch.data(nbchan+6:nbchan+8,:,:);
        %     blink_idx = tar_epoch.data(nbchan+9,:,:);
        %     dataLose_idx = tar_epoch.data(nbchan+10,:,:);
        eyeDirect = tar_epoch.data(nbchan+11:nbchan+13,:,:);
        % calculate GIP distance
        dist_gip2box = dist2Box(GIP, headLoc, tar_lib, tar_direct);

        % calculate rotation
        [headAng, headRot, headAngDiff, headAngCumsum] = cal_rot(headDirect, tar_direct, tar_epoch.srate);
        [eyeAng, eyeRot, eyeAngDiff, eyeAngCumsum] = cal_rot(eyeDirect, tar_direct, tar_epoch.srate);
        [gipAng, gipRot, gipAngDiff, gipAngCumsum] = cal_rot(GIP, tar_direct, tar_epoch.srate);

        % save data to lib
        dist_lib{e_i,subj_i,cond_i} = dist_gip2box;
        ang_lib{e_i,subj_i,cond_i} = {headAng,eyeAng,gipAng};
        angDiff_lib{e_i,subj_i,cond_i} = {headAngDiff,eyeAngDiff,gipAngDiff};
    end
end

disp('Done')
com = 'Dimension: [std_epoch,dev_epoch,gip_std,gip_dev] * subject * condition. Ang_lib cell: [head, eye, gip]';

%%
for j = 1:21
    current_data = data(:, (j-1)*40000+1:j*40000+1);
    Fs = round(1/mean(diff(stream.time_stamps)));
    L = length(current_data);   
    f = Fs*(0:(L/2))/L;

    %%
    Y = fft(current_data);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    %% average with 1Hz bin
    % ignore frequecy higher than 20Hz
    f_right = 40;
    plt_f = f(f<=f_right);
    plt_data = P1(1:length(plt_f));
    % 1Hz average
    compress_rate = 1;
    compress_f = compress_rate:compress_rate:f_right;
    bin = discretize(plt_f,[compress_f-compress_rate/2, compress_f(end)+compress_rate/2]);
    compress_data = zeros(1,length(compress_f));
    for i = 1:length(compress_f)
        compress_data(i) = mean(plt_data(bin==i));
    end

    figure
    plot(compress_f, compress_data,'linewidth',3,'DisplayName',sprintf('%dHz',tarHz));
    grid on
    hold on
    [max_val, max_idx] = max(compress_data);
    plot(compress_f(max_idx),max_val,'ob','DisplayName',sprintf('%dHz',compress_f(max_idx)),'markersize',15);
    title(num2str(j));
    legend
    saveas(gcf, sprintf('%d.png', j))
end


%% TRCA
% tarCh = {test_epoch.up.chanlocs([16:27, 52:59]).labels};
tarCh = {'Oz','O1','O2'};
output = reshape2TRCA(epoch_struct.ring1.stim, tarCh);
[nb_targ,nb_chan,nb_sample,nb_block] = size(output);
fs = test_epoch.up.srate;
num_fbs = 1;
is_ensemble = 0;
accs = zeros(nb_block,1);
labels = 8:11;

% Estimate classification performance
for loocv_i = 1:1:nb_block
    % Training stage 
    traindata = output;
    traindata(:, :, :, loocv_i) = [];
    model = train_trca(traindata, fs, num_fbs);
    
    % Test stage
    testdata = squeeze(output(:, :, :, loocv_i));
    estimated = test_trca(testdata, model, is_ensemble);
    
    % Evaluation 
    is_correct = (estimated==labels);
    accs(loocv_i) = mean(is_correct)*100;
    fprintf('Trial %d: Accuracy = %2.2f%%\n',loocv_i, accs(loocv_i));
end % loocv_i


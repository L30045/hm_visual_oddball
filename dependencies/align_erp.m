function [align_epoch, amp_lib, lat_lib, idx_rm_epoch, reg_mat] = align_erp(epoch, times, varargin)
%% This function performs a multiple linear regression to align ERP.
% principal components calculated from time peirod of interest are used as regressors.
% Input:
%       epoch:          epoch to be aligned (single channel by time points)

%       times:          time points (should be the same length with epoch)

%       t_interest:     time period of interest. Multiple time periods of 
%                       interest should put into a matrix row by row. (default: 0 - 200 ms)

%       ref_epoch:      epoch to calculated regressors (use epoch if reg_epoch is empty)

%       peak_polar:     positive or negative peak to look for. Length should be
%                       the same as the number of t_interest. (default: all positive)
%                       TODO: Use auto correlation

%       nb_PC:          number of principal component to use as regressors (default: 3)

%       mvavg_winlen:   number of epoch in a window when performing moving average
%                       before calculating PCA. (default: 10)

%       t_align:        time period of interest to align (default: mean peak latency within the
%                       first time period of interest.)

%       peak_winlen:    time window used during peak detection. All value
%                       in this time window should be smaller than peak value. (default:20 ms)
%
% Output:
%       align_epoch:    aligned epoch
%       amp_lib:        peak amplitude of each epoch (trial by number of t_interest)
%       lat_lib:        peak latency of each epoch (trial by number of t_interest)
%       idx_rm_epoch:   index of removed epochs
%       reg_mat:        regressors matrix (times by number of regressors)


%% parameter setting
p = inputParser;
p.KeepUnmatched = true;
addRequired(p,'epoch');
addRequired(p,'times');
addOptional(p,'t_interest',[0 200]) % ms (time period of interest)
addOptional(p,'ref_epoch',[])
addOptional(p,'peak_polar',[]) % positve/ negative peak to look for
addOptional(p,'nb_PC',3) % number of principal component to use as regressors
addOptional(p,'mvavg_winlen',10) % number of epoch in a window when performing moving average
%       before calculating PCA.
addOptional(p,'t_align',[]) % time point to align
addOptional(p,'peak_winlen',20) % ms (time window used during peak detection. All value
%       in this time window should be smaller than peak value.)
parse(p,epoch,times,varargin{:})

epoch = p.Results.epoch;
times = p.Results.times;
t_interest = p.Results.t_interest;
ref_epoch = p.Results.ref_epoch;
peak_polar = p.Results.peak_polar;
nb_PC = p.Results.nb_PC;
mvavg_winlen = p.Results.mvavg_winlen;
t_align = p.Results.t_align;
peak_winlen = p.Results.peak_winlen;
% transpose epoch matrix if it is not trial by time
if size(epoch,2)~=length(times)
    epoch = epoch';
end
% assign reg_epoch
if isempty(ref_epoch)
    ref_epoch = epoch;
end
% assign peak_polar
if isempty(peak_polar)
    peak_polar = true(size(t_interest,1),1);
end
    
    
%% Finding PC regressors
% perform moving average on single epoch matrix to reduce noises for the
% following PCA process
random_group = randperm(size(ref_epoch,1));
sort_epoch = ref_epoch(random_group,:);
mvavg_idx = discretize((1:size(sort_epoch,1))-0.5,ceil(size(sort_epoch,1)/mvavg_winlen));
mvavg_epoch = zeros(max(mvavg_idx), size(sort_epoch,2));
for i = 1:max(mvavg_idx)
    mvavg_epoch(i,:) = mean(sort_epoch(mvavg_idx==i,:));
end

% segement out data within time period of interest
reg_mat = zeros(length(times),size(t_interest,1)*nb_PC);
for t_i = 1:size(t_interest,1)
    tmp_t = t_interest(t_i,:);
    idx_t = times <= tmp_t(2) & times >= tmp_t(1);
    pca_epoch = mvavg_epoch(:,idx_t);
    % PCA
    [pc_act,~] = pca(pca_epoch',nb_PC);
    % zero paddling
    pc_regressors = zeros(nb_PC, length(idx_t));
    pc_regressors(:,idx_t) = pc_act;
    % regressors matrix
    reg_mat(:,(t_i-1)*nb_PC+1:t_i*nb_PC) = pc_regressors';
end
reg_mat = [ones(length(times),1), reg_mat];

%% Multiple Linear Regression with dispertion term
reg_epoch = zeros(size(epoch));
for i = 1:size(epoch,1)
    single_epoch = epoch(i,:)';
	b = regress(single_epoch, reg_mat);
    reg_epoch(i,:) = reg_mat*b;
end

% extract time period of interest
reg_epoch = reg_epoch(:,idx_t);

%%===== fix here =====
%%===== fix here =====

% calculate the amplitude and latency
amp_lib = zeros(size(reg_epoch,1),size(t_interest,1));
lat_lib = zeros(size(amp_lib));
t_resolution = mean(diff(times));
for t_i = 1:size(t_interest,1)
    if peak_polar(t_i)
        tmp_p = 1;
    else
        tmp_p = -1;
    end
    tmp_t = t_interest(t_i,:);
    idx_t = times <= tmp_t(2) & times >= tmp_t(1);
    plt_t = times(idx_t);
    t_win = round(peak_winlen/2/t_resolution); %(left/right 10 ms)
    for i = 1:size(reg_epoch,1)
        tmp_epoch = reg_epoch(i,:)*tmp_p;
        [amp, lat] = max(tmp_epoch);
        % remove no peak epoch
        if ~(lat==1 || lat==size(reg_epoch,2)) && all(tmp_epoch(1,max([1 lat-t_win]):min([lat+t_win, size(reg_epoch,2)]))<=amp)
            amp_lib(i,t_i) = amp*tmp_p;
            lat_lib(i,t_i) = plt_t(lat);
        else
            amp_lib(i,t_i) = NaN;
            lat_lib(i,t_i) = NaN;
        end
    end
end

%%===== fix here =====
%%===== fix here =====

%% align epochs
rng(1);
if isempty(t_align)
    t_align = 1;
end
t_delay= lat_lib(:,t_align) - round(nanmean(lat_lib(:,t_align)));
% remove epoch with t_delay more than 50ms
lat_lib(abs(t_delay)>50,t_align) = NaN;

align_epoch = zeros(size(epoch));

for i = 1:size(reg_epoch,1)
    if ~isnan(lat_lib(i,t_align))
        if t_delay(i) > 0
            align_epoch(i,1:end-t_delay(i)+1) = epoch(i,t_delay(i):end);
        elseif t_delay(i) == 0
            align_epoch(i,:) = epoch(i,:);
        else
            align_epoch(i,abs(t_delay(i)):end) = epoch(i,1:end+t_delay(i)+1);
        end
    end
end
% remove epoch with no peak in the t_align period
idx_rm_epoch = isnan(t_delay);

align_epoch(idx_rm_epoch,:) = [];

end
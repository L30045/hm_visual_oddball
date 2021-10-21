function output = linear_interp_data(input, int_p, rm_idx)
%% linear interpret the data to a specific time scale
% Input:
%   input: data to interpret
%   int_p: target time scale to project on
%   rm_idx: input points to ignore
% Output:
%   outpu: interpretted data
% e.g.
%   input = [10, 2, 5, 5, 9; 3, 5, 2, 6, 10];
%   int_p = [1,0,0,1,0,0,0,1,0,0];
%   rm_idx = [1, 5];

%%
output = zeros(size(input,1),length(int_p));
% remove points to ignore
input(:,rm_idx) = [];
int_p = find(int_p);
output(:,int_p(1)) = input(:,1);
floor = int_p(1)+1;
for d_i = 2:size(input,2)
    % floor data
    f_data = real(input(:,d_i-1));
    % ceiling data
    c_data = real(input(:,d_i));
    % linear interpretaion step
    step = real((c_data-f_data)./(int_p(d_i)-floor+1));
    for i_r = 1:size(input,1)
        if isnan(f_data(i_r)) || isnan(c_data(i_r))
            % if f_data is missing, interpret c_data at the time point only
            output(:,int_p(d_i)) = c_data(i_r);
        else
            if abs(step(i_r)) < 1e-04 % tolerance small difference
                output(i_r,floor:int_p(d_i)) = c_data(i_r);
            else
                output(i_r,floor:int_p(d_i)) = f_data(i_r)+step(i_r):step(i_r):c_data(i_r);
            end
        end
    end
    floor = int_p(d_i)+1;
end
output(:,floor:end) = repmat(c_data, 1, size(output,2)-floor+1);

end
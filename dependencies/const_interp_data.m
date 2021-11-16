function output = const_interp_data(input, int_p, rm_idx)
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
%   output = [2,2,2,5,5,5,5,5,5,5;
%             5,5,5,2,2,2,2,6,6,6];


%%
output = zeros(size(input,1),length(int_p));
% remove points to ignore
input(:,rm_idx) = [];
int_p = find(int_p);
for d_i = 1:size(input,2)-1
    output(:,int_p(d_i):int_p(d_i+1)-1) = repmat(input(:,d_i), 1, int_p(d_i+1)-int_p(d_i));
end
output(:,int_p(end):end) = repmat(input(:,end), 1, size(output,2)-int_p(end)+1);

end
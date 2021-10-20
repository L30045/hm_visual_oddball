function output = linear_interp_data(input, int_p, rm_idx)
output = zeros(1,length(int_p));
floor = 1;
int_p = find(int_p);
for d_i = 2:size(input,2)
    if ~ismember(d_i,rm_idx)
        f_data = real(input(d_i-1));
        c_data = real(input(d_i));
        step = real((c_data-f_data)/(int_p(1)-floor));
        if isnan(f_data)
            output(int_p(1)) = c_data;
        elseif step==0
            output(floor:int_p(1)) = c_data;
        else
            output(floor:int_p(1)) = f_data:step:c_data;
        end
        floor = int_p(1)+1;
        int_p(1) = [];
    end
end
end
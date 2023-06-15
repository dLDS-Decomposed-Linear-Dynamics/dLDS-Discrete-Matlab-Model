function vargout = create_reco( F_cell, c_mat, x0, cum_or_ada, real_dyn,addi,bias_vec, order, x_store, to_plot )

% F_cell: cell array of the transport operators;
% c_mat: k X time-1 matrix of coefficients
% x0: k vector (1d)
% cum_or_ada: if cum - use the reconstructed values; if ada - use the right value   
% real dyn: k X time dynamics
% addi = 0
if nargin <5, real_dyn =[]; end
if nargin < 4, cum_or_ada = 'cum';end
if nargin < 6 || isempty(addi), addi = 1; end
if length(size(real_dyn)) == 3, real_dyn = permute(real_dyn, [1,3,2]); end
if nargin < 7 || isempty(bias_vec), bias_vec = zeros(size(real_dyn,1), 1); end
if nargin < 8 || isempty(order), order = 1; end
if nargin < 9 || isempty(x_store), x_store= {}; end
if nargin < 10 || isempty(to_plot), to_plot= false; end
x0 = x0(:);
if strcmp(cum_or_ada, 'cum') || isempty(real_dyn)
    x_vec = [];
    x_prev = x0;
    x_vec(:,1) = x0;
    for time_point = 1:size(c_mat,2)-1
        cur_c = c_mat(:,time_point+addi);
        fx = [];
        for f_num = 1:length(F_cell)
            f_cur = F_cell{f_num};
            fx(:,end+1) = f_cur*x_prev;        
        end
        x_next = fx*cur_c;
        x_vec(:,end+1) = x_next(:)+bias_vec(:);
        x_prev = x_next;
    end
else %if order ==1
    x_vec = [];
    x_vec(:,1) = x0;
    x_prev = x0;
    
    for time_point = 1:size(c_mat,2)-1
        cur_c = c_mat(:,time_point+addi);
        fx = [];
        for f_num = 1:length(F_cell)
            f_cur = F_cell{f_num};
            %disp(size(x_prev))
            fx(:,end+1) = f_cur*x_prev;

        end
        x_next = fx*cur_c;      
        x_vec(:,end+1) = x_next(:);
        if time_point <= size(real_dyn,2)-1-addi
            x_prev = real_dyn(:,time_point+addi+1);
        end
    end
if order > 1
    %x_vec = {};
    x_store{order} = x_vec;
    x_store = create_reco( F_cell, c_mat, x0, cum_or_ada, x_vec,addi,bias_vec, order-1, x_store);
end
   



end
if isempty(x_store) && order == 1

    vargout = x_vec;
else
    x_store{order} = x_vec;
    vargout = x_store;
    %vargout = flip(x_store);
end  


end
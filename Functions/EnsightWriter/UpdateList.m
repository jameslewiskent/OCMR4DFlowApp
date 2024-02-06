function [out_list, list_ind] = UpdateList(in_list, var)
% Aaron Hess, University of Oxford
% 2016
% Check if var is part of in_list, if not add it to out_list and return its
% index in list
% list is a string cell list or a number array
    out_list = in_list;

    if(iscell(in_list))
        % deal as string
        var = char(var);  % has to be a char
        if(size(var,1) > 1)  % ensure a string not a list
            var = var';
        end
        if(isempty(var))
            var='na';
        end
        list_ind = find(strcmp(var,in_list));
        if(isempty(list_ind))
            out_list(end+1) = {var};
            list_ind = length(out_list);
        end

    else % deal as array of numbers
        list_ind = find(in_list == var);
        if(isempty(list_ind))
            out_list(end+1) = var;
            list_ind = length(out_list);
        end
    end


end
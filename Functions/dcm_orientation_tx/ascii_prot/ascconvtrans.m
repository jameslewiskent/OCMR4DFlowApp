function prot = ascconvtrans(in_string)
% Aaron Hess
% University of Oxford
% March 2012
% Function to pars a siemens text protocol into a matlab structure

    global fid;
    global remain_text;
    
    fid = fopen(in_string);
    
    if(fid == 0)
        fprintf('\nunable to open file, parzing string instead\n');
        remain_text = in_string;
    end
    
    
    % look for '### ASCCONV BEGIN ###'
    while (~strcmp(GetLine(),'### ASCCONV BEGIN ###')) && (~IsEmpty())

    end

    if(IsEmpty())
        error('File not compatible')
    end

    prot = prot_empty();  

    line = GetLine();

    % for each line untill '### ASCCONV END ###'
    while ((~IsEmpty()) && (~strncmp(line,'### ASCCONV END ###',19)))
       % Read line untill space
       % seperate variable name and value
        [var_name, value] = strtok(line, ' = ');

        value = sscanf(value,' = %s');
        %depth = length(strfind(var_name, '.')); %how many fields deep does this go
        % tokenize based on '.'
        depth = (length(strfind(var_name, '.'))+1);
        fieldnames = cell(1,depth);
        index = cell(1,depth);
        for i = 1:depth
            [a, var_name] = strtok(var_name,'.');
            %tt  = a;
            %seperate name and array index
            [a, n] = sscanf(a,'%[a-z,A-Z][%d]');
            fieldnames{i} = char(a(1:end-n+1)');
            if n == 2
                index{i} = a(end) + 1; 
            else
                index{i} = 1;
            end
        end
        switch (depth)
            case 1
            prot.(fieldnames{1})(index{1}) = {value};
            case 2
            prot.(fieldnames{1})(index{1}).(fieldnames{2})(index{2}) = {value};
            case 3
            prot.(fieldnames{1})(index{1}).(fieldnames{2})(index{2}).(fieldnames{3})(index{3}) =  {value};
            case 4
            prot.(fieldnames{1})(index{1}).(fieldnames{2})(index{2}).(fieldnames{3})(index{3}).(fieldnames{4})(index{4}) =  {value};
            case 5
            prot.(fieldnames{1})(index{1}).(fieldnames{2})(index{2}).(fieldnames{3})(index{3}).(fieldnames{4})(index{4}).(fieldnames{5})(index{5}) =  {value};
            case 6
            prot.(fieldnames{1})(index{1}).(fieldnames{2})(index{2}).(fieldnames{3})(index{3}).(fieldnames{4})(index{4}).(fieldnames{5})(index{5}).(fieldnames{6})(index{6}) =  {value};
        end
        % seperate off arrage index [*] on each token
        % 
        % buld into structures by the names 
        line = GetLine();
    end

    fclose(fid);
end


function line = GetLine()
    global remain_text;
    if(fid==0)
        [line,remain_text] = strtok(remain_text,'\n');
    else
        line = fgetl(fid);
    end
end

function is_end = IsEmpty()
    global fid;
    global remain_text;
    if(fid==0)
        [test_line,~] = strtok(remain_text,'\n');
        is_end = isempty(test_line);
    else
        is_end = feof(fid);
    end
end


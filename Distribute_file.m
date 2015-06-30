function[] = Distribute_file()

% Create new directories
clc;
parent_dir = '/Users/xyao/matlab/code-private/fig_organized_temp/';
original_dir = '/Users/xyao/matlab/code-private/fig_temp/';
not_done = 1;
while(not_done)
    dirname_temp = input('directory name:');
    mkdir(parent_dir, dirname_temp);
    not_done = input('create another directory? Yes:1 No:0\n');
end

    
% distribute files
not_done = 1;
files = dir(original_dir);
while(not_done)    
    if(isempty(files))
        not_done = 0;
    elseif(files(1).isdir)
        files(1) = [];
    else
        filename_temp = files(1).name;
        dirname_temp = dir(parent_dir);
        dirname = cell(0,0);
        for i = 1:length(dirname_temp)
            if dirname_temp(i).name(1)~='.'
                dirname = [dirname dirname_temp(i).name];
            end
        end
        dirnum = [];
        for i = 1:length(dirname)
            dirnum = [dirnum '  ' num2str(i) ':' dirname{i}];
        end
        fprintf([dirnum '\n'])
        cd(original_dir);
        open(filename_temp);
        dirn = input('number of target directory:');
        destination = [parent_dir dirname{dirn}];
        movefile(filename_temp, destination)
        files(1) = [];
    end
end

        
        
    
    

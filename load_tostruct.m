function struct_output = load_tostruct(filename)
struct_output = struct;
load(filename)
Varname = who;
for i = 1:length(Varname)
    if ~strcmp(Varname{i}, 'struct_output')
        struct_output = setfield(struct_output, Varname{i}, eval(Varname{i}));
    end
end

end

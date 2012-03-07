function ret = eventlog_parser(eventlogfile)
f = fopen(eventlogfile);
res = textscan(f, '%s %f %f %f %f %f %f', 'Delimiter', '\t');
fclose(f);
names= res{1};
times = res{2}; 
minimum = res{3};
average = res{4};
maximum = res{5};
total = res{6};
rate = res{7};

uniquenames = unique(names);
numentries = length(names);
ret = {};
if (isempty(uniquenames))
    return
end


for i = 1:length(uniquenames)
    ret{i} = struct('name', [], ...
                'times', [], ...
                'minimum', [], ...
                'average', [], ...
                'maximum', [], ...
                'total', [], ...
                'rate', []);

    ret{i}.name = uniquenames{i};
    for j = 1:numentries
        if (strcmp(uniquenames{i}, names{j}))
            ret{i}.times = [ret{i}.times, times(j)];
            ret{i}.minimum = [ret{i}.minimum, minimum(j)];
            ret{i}.average = [ret{i}.average, average(j)];
            ret{i}.maximum = [ret{i}.maximum, maximum(j)];
            ret{i}.total = [ret{i}.total, total(j)];
            ret{i}.rate = [ret{i}.rate , rate(j)];
        end
    end
end
end
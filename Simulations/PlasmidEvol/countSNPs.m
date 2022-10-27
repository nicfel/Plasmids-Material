clear
fastas = dir('out/*nexus');
f = fopen('snpcounts.txt', 'w');
fprintf(f,'run\tplasmid\tsnpcount\n');

for i = 1 : length(fastas)
    g = fopen(['out/' fastas(i).name]);
    seq = '';
    while ~feof(g)
        line = fgets(g);
        if contains(line, 'matrix')
            while ~feof(g)
                line = fgets(g);
                line = strtrim(line);
                if startsWith(line, 't')
                    tmp = strrep(strsplit(line), ';','');
                    seq = [seq;tmp{2}];
                end
            end
        end        
    end
    snp_count = 0;
    for j = 1 : length(seq)
        l = length(unique(seq(:,j)));
        if l>1
            snp_count=snp_count+1;
        end
    end
    tmp = strsplit(fastas(i).name, '.');
    fprintf(f, '%s\t%s\t%d\n',strrep(tmp{1}, 'sim_',''), tmp{2}, snp_count);
end
fclose(f);
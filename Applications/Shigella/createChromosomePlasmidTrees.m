% builds raxml trees from the alignment for chromosome plasmid pairs
clear
system('rm -r tangletrees')
system('mkdir tangletrees')

plasmids = {'Sson_NC_009345_spA_over70.full','Sson_NC_009345_4AMRgenes.full','Sson_NC_009345_strABsul2.full'};

isolates = cell(0,0);
for i = 1:length(plasmids)
    aln = fastaread(['data/Sonnei/' plasmids{i} '.aln']);
    for j = 1:length(aln)
        aln(j).Sequence = strrep(aln(j).Sequence, 'X', 'N');
    end
    % remove the isolate with name Reference
    for j = 1:length(aln)
        if contains(aln(j).Header, 'Reference')
            aln(j) = [];
        end
    end
    fastawrite(['tangletrees/' plasmids{i} '.aln'], aln);

    % combine the headers into a single cell array
    for j = 1:length(aln)
        isolates{end+1} = aln(j).Header;
    end
end


reference_isolates = unique(isolates);
% make a chromsome alignment using the reference isolates
aln = fastaread('data/Sonnei/Sson_NC_007384_18122023_cleanGubbinsV241.filtered_polymorphic_sites.snpsites.fasta');
% get the isolate names
clear chromsome
c=1;
for i = 1:length(aln)
    if any(strcmp(aln(i).Header, reference_isolates))
        chromsome(c) = aln(i);
        c=c+1;
    end
end
% write the chromosome alignment to a file
fastawrite('tangletrees/chromosome.aln', chromsome);

% build a raxml-ng tree for each alignment in the tangletrees directory
files = dir('tangletrees/*.aln');
for i = 1:length(files)
    % check if the final tree already exists
    if exist(['tangletrees/' files(i).name(1:end-4) '_divergence_tree.nexus'], 'file') == 2
        continue
    end
    system(['/opt/homebrew/bin/iqtree2 -nt 11 -s tangletrees/' files(i).name ' -m GTR --prefix tangletrees/' files(i).name(1:end-4)]);
    % run treetime on the tree and alignment
    if contains(files(i).name, 'chromosome')
        system(['/opt/homebrew/Caskroom/miniconda/base/envs/plasmids/bin/python -m treetime --tree tangletrees/' files(i).name(1:end-4)...
            '.treefile --aln tangletrees/' files(i).name...
            ' --dates data/Sonnei_MonthYear_07082022.csv'...
            ' --outdir tangletrees/ ']);
        % rename the file divergence_tree.nexus to the name of the alignment
    else
        system(['/opt/homebrew/Caskroom/miniconda/base/envs/plasmids/bin/python -m treetime --tree tangletrees/' files(i).name(1:end-4)...
            '.treefile --aln tangletrees/' files(i).name...
            ' --dates data/Sonnei_MonthYear_07082022.csv'...
            ' --clock-rate 0.000003 --outdir tangletrees/ ']);
    end
    system(['mv tangletrees/divergence_tree.nexus tangletrees/' files(i).name(1:end-4) '_divergence_tree.nexus']);
end





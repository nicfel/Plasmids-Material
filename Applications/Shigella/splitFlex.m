flex = fastaread('data/Aln_Trees/Flexneri_ref/Sflex_NC_004337_26082022clean_gubbinsv241.filtered_polymorphic_sites.fasta');
pdist = zeros(length(flex), length(flex));
for a = 1 : length(flex)
    disp(a)
    for b = a+1:length(flex)
        indices = find(flex(a).Sequence~=flex(b).Sequence);
        red_indices = find(flex(a).Sequence(indices)=='-' | flex(b).Sequence(indices)=='-');
        gaps = find(flex(a).Sequence=='-' | flex(b).Sequence=='-');

        pdist(a,b) = (length(indices)-length(red_indices))/(length(flex(a).Sequence)-length(gaps));
        pdist(b,a) = pdist(a,b);
    end
end
members = cell(0,0);
already_clustered = [];
min_dist=0.3;
c = 1;
cl_size = 0;
for i = 1:length(flex)
    if ~ismember(i, already_clustered)
        newmembers = find(pdist(i, :)<min_dist);
        members{c} = newmembers;
        if sum(ismember(newmembers, already_clustered))>0
            error('alal');
        end
        while ~isempty(newmembers)
            tmp = [];
            for j = 1 : length(newmembers)
                tmp = [tmp,find(pdist(newmembers(j), :)<min_dist)];
            end
            tmp = unique(tmp);
            newmembers = tmp(~ismember(tmp,members{c}));

            if sum(ismember(newmembers, already_clustered))>0
                error('alal');
            end
            members{c} = [members{c} newmembers];
        end
        cl_size = cl_size+length(members{c});
        already_clustered = [already_clustered, members{c}];
        c=c+1;    
    end
end   

fasta1 = flex(members{1});
fasta2 = flex(members{2});    

% get all constant sites from fasta1 and fasta2
seq='';
for i=1:length(fasta1)
   seq=[seq;fasta1(i).Sequence];
end
const1 = false(1,length(fasta1(1).Sequence));
for i=1:length(seq)
    u = unique(seq(:,i));
    u(u=='N') = [];
    u(u=='-') = [];
    if length(u)<2
        const1(i) = true;
    end
end
seq='';
for i=1:length(fasta2)
   seq=[seq;fasta2(i).Sequence];
end

const2 = false(1,length(fasta2(1).Sequence));
for i=1:length(seq)
    u = unique(seq(:,i));
    u(u=='N') = [];
    u(u=='-') = [];
    if length(u)<2
        const2(i) = true;
    end
end

for i = 1: length(fasta1)
    fasta1(i).Sequence = fasta1(i).Sequence(~const1);
end
for i = 1: length(fasta2)
    fasta2(i).Sequence = fasta2(i).Sequence(~const2);
end


delete('data/Aln_Trees/FlexSon_Ss046/Flex1_NC_007384.aln');
fastawrite('data/Aln_Trees/Flexneri_ref/Flex1_NC_007384.aln',fasta1);

delete('data/Aln_Trees/FlexSon_Ss046/Flex2_NC_007384.aln');
fastawrite('data/Aln_Trees/Flexneri_ref/Flex2_NC_007384.aln',fasta2);



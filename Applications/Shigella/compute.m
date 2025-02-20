NC_007384 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/Sson_NC_007384_29042022clean_gubbinsv241.filtered_polymorphic_sites.fasta');
NC_007385_sub = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_007385_plasmid_SonneiOnly_Over70_0408202.aln');
NC_007385 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_007385_plasmid_SonneiOnly_Over50_2408202_full_clean30092022.aln');

NC_009345 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009345_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln');
NC_009346 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009346_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln');
NC_009347 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009347_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln');



clear str perc
for i = 1 : length(NC_007385)
    str(i,:) = strrep(NC_007385(i).Sequence, 'N','-');
    perc(i) = 1-sum(str(i,:)=='-')/length(str(i,:));    
end
str_new = str(perc>0.7,:);


f = fopen('shigella_lengthrate.tsv', 'w');
fprintf(f, 'plasmid\tsnps\tlengthforrate\n');
fprintf(f, 'NC_007384\t%d\t%d\n',length(NC_007384(1).Sequence), 4825265);
fprintf(f, 'NC_007385\t%d\t%d\n', length(NC_007385_sub(1).Sequence), sum(sum(str_new=='-')==0));
fprintf(f, 'NC_009345\t%d\t%d\n', length(NC_009345(1).Sequence),length(NC_009345(1).Sequence));
fprintf(f, 'NC_009346\t%d\t%d\n', length(NC_009345(1).Sequence),length(NC_009346(1).Sequence));
fprintf(f, 'NC_009347\t%d\t%d\n', length(NC_009345(1).Sequence),length(NC_009347(1).Sequence));


LN624486_sub = fastaread('data/Aln_Trees/LN624486_MDR_plasmid_aln/LN624486_SonFlex_Over70_26082022.aln');
LN624486 = fastaread('data/Aln_Trees/LN624486_MDR_plasmid_aln/LN624486_SonFlex_Over50_26082022_full_clean30092022.aln');

clear str perc
for i = 1 : length(LN624486)
    str(i,:) = strrep(LN624486(i).Sequence, 'N','-');
    perc(i) = 1-sum(str(i,:)=='-')/length(str(i,:));    
end
str_new = str(perc>0.7,:);

fprintf(f, 'LN624486\t%d\t%d\n',length(LN624486_sub(1).Sequence),  sum(sum(str_new=='-')==0));

length(LN624486_sub(1).Sequence)/sum(sum(str_new=='-')==0);


fclose(f);
clear
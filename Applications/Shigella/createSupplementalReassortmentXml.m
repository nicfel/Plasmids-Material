rng(187989)

clear

weights_name = {'equal','weighted'};


%% Shigella sonnei
meta = importdata('data/Sonnei_MonthYear_07082022.txt');

names = {
        'Sson_NC_009345_spA_over70.full.aln';
        'Sson_NC_009345_strABsul2.full.aln';
        'Sson_NC_009345_4AMRgenes.full.aln';
        'Sson_NC_007384_18122023_cleanGubbinsV241.filtered_polymorphic_sites.snpsites.fasta'};




fnames = {'spA', 'strABsul2', '4AMRgenes', 'chromosome'};
spa = struct2table(fastaread(['data/Sonnei/' names{1}]));
str = struct2table(fastaread(['data/Sonnei/' names{2}]));
amr = struct2table(fastaread(['data/Sonnei/' names{3}]));
cr = struct2table(fastaread(['data/Sonnei/' names{3}]));

ids = meta.textdata(2:end,1);
times = strrep(meta.textdata(2:end,2), '_','-');


has_plasmid = -1*ones(length(ids), 6);
for i = 1 : length(ids)
    if ismember(ids{i}, spa.Header)
        has_plasmid(i,1) = find(ismember(spa.Header,ids{i}));
    end
    if ismember(ids{i}, str.Header)
        has_plasmid(i,2) = find(ismember(str.Header,ids{i}));
    end
    if ismember(ids{i}, amr.Header)
        has_plasmid(i,3) = find(ismember(amr.Header,ids{i}));
    end
    if ismember(ids{i}, cr.Header)
        has_plasmid(i,4) = find(ismember(cr.Header,ids{i}));
    end
end
% das

weights = (has_plasmid(:,1)>0) + (has_plasmid(:,2)>0) + (has_plasmid(:,3)>0);

segs_files(1).name = names{1};
segs_files(2).name = names{2};
segs_files(3).name = names{3};
segs_files(4).name = names{4};

use_seqs = sort(find(weights==3));

for i = 1:4
    fasta = fastaread(['data/Sonnei/' segs_files(i).name]);
    clear Data
    for j = 1 : length(use_seqs)
        ind = has_plasmid(use_seqs(j),i);
        h = fasta(ind).Header;
        header ='';
        for k = 1: length(meta.textdata) 
            if strcmp(h, meta.textdata{k})
                header = [meta.textdata{k,1} '|' strrep(meta.textdata{k,2}, '_','-') '-01'];
            end
        end
        Data(j).Header = header;
        Data(j).Sequence =strrep(fasta(ind).Sequence, 'X','N');
    end
    delete(['spa/',fnames{i} '.fasta'])
    fastawrite(['spa/',fnames{i} '.fasta'],  Data);
end


rng(187989)

clear


%% Shigella, create 
meta = importdata('data/Sonnei_MonthYear_07082022.txt');

names = {'Sson_NC_007384_18122023_cleanGubbinsV241.filtered_polymorphic_sites.snpsites.fasta';
        'Sson_NC_007385_VP_over60.filtered95.aln';
        'Sson_NC_009345_spA_over70.full.aln';
        'Sson_NC_009346_spB_over90.full.aln';
        'Sson_NC_009347_spC_over90.full.aln'};

NC_007384 = struct2table(fastaread(['data/Sonnei/' names{1}]));
NC_007385 = struct2table(fastaread(['data/Sonnei/' names{2}]));
NC_009345 = struct2table(fastaread(['data/Sonnei/' names{3}]));
NC_009346 = struct2table(fastaread(['data/Sonnei/' names{4}]));
NC_009347 = struct2table(fastaread(['data/Sonnei/' names{5}]));

ids = meta.textdata(2:end,1);
times = strrep(meta.textdata(2:end,2), '_','-');


has_plasmid = -1*ones(length(ids), 5);
for i = 1 : length(ids)
    if ismember(ids{i}, NC_007384.Header)
        has_plasmid(i,1) = find(ismember(NC_007384.Header,ids{i}));
    end
    if ismember(ids{i}, NC_007385.Header)
        has_plasmid(i,2) = find(ismember(NC_007385.Header,ids{i}));
    end
    if ismember(ids{i}, NC_009345.Header)
        has_plasmid(i,3) = find(ismember(NC_009345.Header,ids{i}));
    end
    if ismember(ids{i}, NC_009346.Header)
        has_plasmid(i,4) = find(ismember(NC_009346.Header,ids{i}));
    end
    if ismember(ids{i}, NC_009347.Header)
        has_plasmid(i,5) = find(ismember(NC_009347.Header,ids{i}));
    end     
end
% das

weights = (has_plasmid(:,2)>0) + (has_plasmid(:,3)>0) + (has_plasmid(:,4)>0) + (has_plasmid(:,5)>0);

weights(:) = 1;            


indices = [];
for i = 1 : 400
    use_weights = weights;
    use_weights(indices) = 0;
    indices = [indices, randsample(length(weights), 1, true, use_weights)];
end
use_seqs = sort(indices);


segs_files(1).name = names{1};
segs_files(2).name = names{2};
segs_files(3).name = names{3};
segs_files(4).name = names{4};
segs_files(5).name = names{5};

i=1;
fasta = fastaread(['data/Sonnei/' segs_files(i).name]);
clear Data

true_plasmids = [];

for j = 1 : length(use_seqs)
    ind = has_plasmid(use_seqs(j),i);
    h = fasta(ind).Header;
    header ='';
    for k = 1:length(meta.textdata) 
        if strcmp(h, meta.textdata{k})
            header = [meta.textdata{k,1} '|' strrep(meta.textdata{k,2}, '_','-') '-01'];
        end
    end
    true_plasmids(j,:) = has_plasmid(use_seqs(j),2:end)>-1;

    Data(j).Header = header;
    Data(j).Sequence =strrep(fasta(ind).Sequence, 'X','N');
end

rand_ind = randsample(size(true_plasmids,1),size(true_plasmids,1));
random_tip = true_plasmids(rand_ind,:);
for j = 1 : length(Data)
    Data(j).Header = [Data(j).Header sprintf('|%d',true_plasmids(j,:)) sprintf('|%d',random_tip(j,:))];
end

delete(['supplementalxmls/dta.fasta']);
fastawrite(['supplementalxmls/dta.fasta'],  Data);

%% fill in template 
f = fopen('supplementalxmls/dta_template.xml');
g = fopen('supplementalxmls/dta.xml', 'w');

names = {'pInv', 'spA', 'spB', 'spC', 'pInv_rand', 'spA_rand', 'spB_rand', 'spC_rand'};

while ~feof(f)
    line = fgets(f);
    if contains(line, 'insert_tips')
        for j = 1 : length(Data)
            Data(j).Header = [Data(j).Header sprintf('|%d',true_plasmids(j,:)) sprintf('|%d',random_tip(j,:))];
            fprintf(g, '\t\t<sequence id="seq_%s" spec="Sequence" taxon="%s" totalcount="4" value="%s"/>\n',Data(j).Header,Data(j).Header,Data(j).Sequence)
        end
    elseif contains(line, 'insert_dates')
        tmp = '';
        for j = 1 : length(Data)
            tmp2 = strsplit(Data(j).Header, '|');
            tmp = [tmp ',' Data(j).Header '=' tmp2{2}];
        end
        fprintf(g, strrep(line, 'insert_dates', tmp(2:end)));
    elseif contains(line, 'insertParams')
        for j = 1:length(names)
            fprintf(g, '\t\t\t\t<parameter id="relativeGeoRates.s:%s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>\n', names{j});
            fprintf(g, '\t\t\t\t<parameter id="traitClockRate.c:%s" spec="parameter.RealParameter" name="stateNode">1.0</parameter>\n', names{j});
        end
    elseif contains(line, 'insert_prior')
                for j = 1:length(names)
                 fprintf(g, '\t\t\t\t<prior name="distribution" x="@relativeGeoRates.s:%s">\n', names{j});
                 fprintf(g, '\t\t\t\t<Gamma name="distr">\n');
                 fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" estimate="false" name="alpha">1.0</parameter>\n');
                 fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" estimate="false" name="beta">1.0</parameter>\n');
                 fprintf(g, '\t\t\t\t</Gamma>\n');
                 fprintf(g, '\t\t\t\t</prior>\n');
                  fprintf(g, '\t\t\t\t<prior name="distribution" x="@traitClockRate.c:%s">\n', names{j});
                 fprintf(g, '\t\t\t\t<LogNormal name="distr" M="1" S="1.25"/>\n');
                 fprintf(g, '\t\t\t\t</prior>\n');

                end
    elseif contains(line, 'insert_likelihood')
        for j = 1:length(names)
            fprintf(g, '\t\t\t\t<distribution id="traitedtreeLikelihood.%s" spec="beastclassic.evolution.likelihood.AncestralStateTreeLikelihood" tag="%s" tree="@Tree.t:dta">\n', names{j}, names{j});
            fprintf(g, '\t\t\t\t<data id="%s" spec="beastclassic.evolution.alignment.AlignmentFromTrait">\n', names{j});
            fprintf(g, '\t\t\t\t<traitSet id="traitSet.%s" spec="beast.base.evolution.tree.TraitSet" taxa="@TaxonSet.dta" traitname="discrete">\n', names{j});
            tmp = '';
            for k = 1:length(Data)
                tmp2 = strsplit(Data(k).Header, '|');
                tmp = [tmp ',' Data(k).Header '=' tmp2{2+j}];
            end
            fprintf(g, '\t\t\t\t%s</traitSet>\n', tmp(2:end));

                     fprintf(g, '\t\t\t\t<userDataType id="traitDataType.%s" spec="beast.base.evolution.datatype.UserDataType" codeMap="0=0,1=1,? = 0 1 " codelength="-1" states="2"/>\n', names{j});
                     fprintf(g, '\t\t\t\t</data>\n');
                     fprintf(g, '\t\t\t\t<siteModel id="geoSiteModel.s:%s" spec="SiteModel" gammaCategoryCount="1">\n', names{j});
                         fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>\n', names{j});
                         fprintf(g, '\t\t\t\t<parameter id="gammaShape.s:%s" spec="parameter.RealParameter" estimate="false" name="shape">1.0</parameter>\n', names{j});
                         fprintf(g, '\t\t\t\t<parameter id="proportionInvariant.s:%s" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n', names{j});
                         fprintf(g, '\t\t\t\t<substModel id="svs.s:%s" spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel" rates="@relativeGeoRates.s:%s" symmetric="false">\n', names{j}, names{j});
                             fprintf(g, '\t\t\t\t<rateIndicator id="rateIndicator.s:%s" spec="parameter.BooleanParameter" estimate="false">true</rateIndicator>\n', names{j});
                             fprintf(g, '\t\t\t\t<frequencies id="traitfreqs.s:%s" spec="Frequencies">\n', names{j});
                                 fprintf(g, '\t\t\t\t<parameter id="traitfrequencies.s:%s" spec="parameter.RealParameter" dimension="2" name="frequencies">0.5</parameter>\n', names{j});
                             fprintf(g, '\t\t\t\t</frequencies>\n');
                         fprintf(g, '\t\t\t\t</substModel>\n');
                     fprintf(g, '\t\t\t\t</siteModel>\n');
                     fprintf(g, '\t\t\t\t<branchRateModel id="StrictClockModel.c:%s" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@traitClockRate.c:%s"/>\n', names{j}, names{j});
                fprintf(g, '\t\t\t\t</distribution>\n');
        end
    elseif contains(line, 'insert_scaler')
        for j = 1:length(names)

            fprintf(g, '\t\t\t\t<operator id="georateScaler.s:%s" spec="ScaleOperator" parameter="@relativeGeoRates.s:%s" scaleAllIndependently="true" scaleFactor="0.99" weight="30.0"/>\n', names{j}, names{j});
            fprintf(g, '\t\t\t\t<operator id="geoMuScaler.c:%s" spec="ScaleOperator" parameter="@traitClockRate.c:%s" scaleFactor="0.9" weight="3.0"/>\n', names{j}, names{j});
        end
    elseif contains(line, 'insert_log')
        for j = 1:length(names)
             fprintf(g, '\t\t\t\t<log idref="relativeGeoRates.s:%s"/>\n', names{j});
             fprintf(g, '\t\t\t\t<log idref="traitClockRate.c:%s"/>\n', names{j});
             fprintf(g, '\t\t\t\t<log id="geoSubstModelLogger.s:%s" spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModelLogger" dataType="@traitDataType.%s" model="@svs.s:%s"/>\n', names{j}, names{j}, names{j});
        end
    elseif contains(line, 'insert_tree')

        for j = 1:length(names)
             fprintf(g, '\t\t\t\t<metadata idref="traitedtreeLikelihood.s:%s"/>\n', names{j});
        end
    else
        fprintf(g, line);
    end





 
end


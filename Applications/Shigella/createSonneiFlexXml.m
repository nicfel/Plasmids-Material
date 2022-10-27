rng(187989)

clear

clockrates = [3.444E-4 8.5E-5 1.662E-5];
for cutoff = 50

    %% Shigella sonnei
    meta1 = importdata('data/Sonnei_MonthYear_07082022.txt');
    meta2 = importdata('data/ShigFlex_01092022_YearMonth.txt');
    
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
    
    % check that there are two groups exactly
    if length(members)~=2
        error('wrong group number');
    end
    fasta1 = flex(members{1});
    fasta2 = flex(members{2});    
    fasta3 = fastaread('data/Aln_Trees/Ss046_Sonnei_ref/Sson_NC_007384_29042022clean_gubbinsv241.filtered_polymorphic_sites.fasta');
    
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
    
    c=1;    
    clear fasta_comb;    
    % combine fastas
    for i = 1 : length(fasta1)
        fasta_comb(c).Header = fasta1(i).Header;
        fasta_comb(c).Sequence = [fasta1(i).Sequence(~const1), repmat('N', 1,sum(~const2)+length(fasta3(1).Sequence))];
        c=c+1;
    end
    for i = 1 : length(fasta2)
        fasta_comb(c).Header = fasta2(i).Header;
        fasta_comb(c).Sequence = [repmat('N', 1,sum(~const1)) fasta2(i).Sequence(~const2) repmat('N',1,length(fasta3(1).Sequence))];
        c=c+1;
    end    
    for i = 1 : length(fasta3)
        fasta_comb(c).Header = fasta3(i).Header;
        fasta_comb(c).Sequence = [repmat('N', 1,sum(~const1)+sum(~const2)) fasta3(i).Sequence];
        c=c+1;
    end
    
    
    
    
    delete('data/Aln_Trees/FlexSon_Ss046/SonFlex_NC_007384_concat.aln');
    fastawrite('data/Aln_Trees/FlexSon_Ss046/SonFlex_NC_007384_concat.aln',fasta_comb);

    NC_007384 = struct2table(fastaread('data/Aln_Trees/FlexSon_Ss046/SonFlex_NC_007384_concat.aln'));
%     NC_007385 = struct2table(fastaread('data/Aln_Trees/FlexSon_Ss046/SonFlex_NC_007385_Over50_23082022.aln'));
    LN624486 = struct2table(fastaread('data/Aln_Trees/LN624486_MDR_plasmid_aln/LN624486_SonFlex_Over70_26082022.aln'));
    
    ids = [meta1.textdata(2:end,1);meta2.textdata(2:end,1)];
    times = strrep([meta1.textdata(2:end,2);meta2.textdata(2:end,2)], '_','-');

    has_plasmid = -1*ones(length(ids), 3);
    for i = 1 : length(ids)
        if ismember(ids{i}, NC_007384.Header)
            has_plasmid(i,1) = find(ismember(NC_007384.Header,ids{i}));
            if NC_007384.Sequence{has_plasmid(i,1)}(1)~='N'
                group(i) = 1;
            elseif NC_007384.Sequence{has_plasmid(i,1)}(end)~='N'
                group(i) = 3;
            else
                group(i) = 2;
            end
        end
        if ismember(ids{i}, LN624486.Header)
            has_plasmid(i,2) = find(ismember(LN624486.Header,ids{i}));
        end
%         if ismember(ids{i}, NC_007385.Header)
%             has_plasmid(i,3) = find(ismember(NC_007385.Header,ids{i}));
%         end

    end

    dont_use_seqs = ~(has_plasmid(:,2)>0);
    
    plasmid_group = group;
    plasmid_group(dont_use_seqs) =-1;

    for i = 1 : max(plasmid_group)
        weights(plasmid_group==i) = 1/sum(plasmid_group==i);
    end
    weights(dont_use_seqs) = 0;
    
    indices = [];
    for i = 1 : 300
        use_weights = weights;
        use_weights(indices) = 0;
        indices = [indices, randsample(length(weights), 1, true, use_weights)];
    end
    use_seqs = sort(indices);

    segs_files(1).name = 'SonFlex_NC_007384_concat.aln';
    segs_files(2).name = 'LN624486_SonFlex_Over70_26082022.aln';
%     segs_files(3).name = 'SonFlex_NC_007385_Over70_23082022.aln';

    segments = {'NC_007384','LN624486'};

%     dsa

    %%
    for r = 0:2
        f = fopen('template.xml');
        g = fopen(['xmls/SonFlex_rep' num2str(r) '.xml'], 'w');

        est_tip_time = cell(0,0);
        while ~feof(f)
            line = fgets(f);
            if contains(line, 'insert_alignment')
                for i = 1 : length(segments)
                    fprintf(g, '\t<data id="%s">\n',segments{i});
                    if startsWith(segs_files(i).name,'LN')
                        fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(i).name]);
                    else
                        fasta = fastaread(['data/Aln_Trees/FlexSon_Ss046/' segs_files(i).name]);
                    end

                    seq_length(i) = length(fasta(1).Sequence);

                    for j = 1 : length(use_seqs)
                        ind = has_plasmid(use_seqs(j),i);
                        if ind==-1
                        else
                            fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                 segments{i}, fasta(ind).Header,...
                                 fasta(ind).Header, fasta(ind).Sequence);
                        end
                    end
                    fprintf(g, '\t</data>\n');
                    if i==1
                        % add dummy alignment for inititalization
                        fprintf(g, '\t<data id="dummy">\n');
                        for j = 1 : length(use_seqs)
                            length_vals = min(sum(~const1),length(fasta3(1).Sequence));
                            start_vals = sum(~const1)+sum(~const2)+1;

                            ind = has_plasmid(use_seqs(j),i);
                            if ind==-1
                            else
                                if group(use_seqs(j))==1
                                    seq = fasta(ind).Sequence(1:length_vals);
                                else
                                    seq = fasta(ind).Sequence(start_vals:(start_vals+length_vals-1));
                                end
                                fprintf(g, '\t\t<sequence id="dummy.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                     fasta(ind).Header,...
                                     fasta(ind).Header, seq);
                            end
                        end
                        fprintf(g, '\t</data>\n');

                    end


                end
            elseif contains(line, 'insert_run_header')
                fprintf(g, '\t\t<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" chainLength="10000000" deltaTemperature="0.00001" heatLikelihoodOnly="true" storeEvery="1000000" chains="2" resampleEvery="10000">\n');            
%                  fprintf(g, '\t\t<run id="mcmc" spec="beast.core.MCMC" chainLength="10000000" storeEvery="1000000">\n');            
            elseif contains(line, 'insert_nr_plasmids')
                fprintf(g, strrep(line, 'insert_nr_plasmids', num2str(length(segments)-1)));

            elseif contains(line, 'insert_taxa')
                fprintf(g, '\t\t<taxonset id="taxonSet" spec="TaxonSet">\n');
                for j = 1 : length(use_seqs)
                    fprintf(g, '\t\t\t<taxon spec="Taxon" id="%s"/>\n', ids{use_seqs(j)});
                end
                fprintf(g, '\t\t</taxonset>\n');
            elseif contains(line, 'insert_taxonsets')
                for i = 2 : length(segments)
                    fprintf(g, '\t\t\t\t<plasmidTaxonSet idref="TaxonSet.%s" spec="TaxonSet"/>\n', segments{i});
                end
            elseif contains(line, 'insert_sampling_times')
                for j = 1 : length(use_seqs)
                    if j < length(use_seqs)
                        fprintf(g, '%s=%s-01,\n', ids{use_seqs(j)}, times{use_seqs(j)});
                    else
                        fprintf(g, '%s=%s-01\n', ids{use_seqs(j)}, times{use_seqs(j)});
                    end
                end  
            elseif contains(line, 'insert_nr_segments')
                fprintf(g, strrep(line, 'insert_nr_segments', num2str(length(segments))));
            elseif contains(line, 'insert_parameters')

                for s = 1 : length(segments)
                    if s==1 
%                         fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" name="stateNode" lower="2.5e-04" upper="4e-04">3.5e-04</parameter>\n',segments{s});

                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
                                if i==3
                                    fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s.%d" name="stateNode" lower="2.5e-04" upper="4e-04">%f</parameter>\n',segments{s}, i, clockrates(1));
                                else
                                    fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s.%d" name="stateNode" lower="2.e-05" upper="1e-04">%f</parameter>\n',segments{s}, i, clockrates(2));
                                end
                                fprintf(g, '\t\t\t\t<parameter id="kappa.s:%s.%d" lower="0.0" name="stateNode">%f</parameter>\n',segments{s},i, lognrnd(0,0.5,1));
                                fprintf(g, '\t\t\t\t<parameter id="gammaShape.s:%s.%d" name="stateNode">%f</parameter>\n',segments{s},i, lognrnd(0,0.5,1));
                                fprintf(g, '\t\t\t\t<parameter id="freqParameter.s:%s.%d" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',segments{s},i);

                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" name="stateNode">0.005</parameter>\n',segments{s});
                        fprintf(g, '\t\t\t\t<parameter id="kappa.s:%s" lower="0.0" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                        fprintf(g, '\t\t\t\t<parameter id="gammaShape.s:%s" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                        fprintf(g, '\t\t\t\t<parameter id="freqParameter.s:%s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',segments{s});
                    end
                end
            elseif contains(line, 'insert_corename')
                tmpline = strrep(line, 'insert_corename.tree',[segments{1} '.tree']);
                fprintf(g, strrep(tmpline, 'insert_corename','dummy'));
            elseif contains(line, 'insert_priors')
                % insert sampling time priors
                for s = 1 : length(segments)
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
%                                 fprintf(g, '\t\t\t\t\t\t\t\t<prior id="MutPrior.s:%s.%d" name="distribution" x="@mutationRate.s:%s.%d">\n', segments{s}, i, segments{s}, i);
%                                 fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistribution.%s.%d"  meanInRealSpace="true" M="0.0005" S="2" name="distr"/>\n', segments{s}, i);
%                                 fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                                fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s.%d" name="distribution" x="@kappa.s:%s.%d">\n', segments{s},i, segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1.%d" name="distr" M="1.0" S="1.25"/>\n', segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                                fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s.%d" name="distribution" x="@gammaShape.s:%s.%d">\n', segments{s},i, segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s.%d" name="distr"/>\n', segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');

                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s" name="distribution" x="@kappa.s:%s">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1" name="distr" M="1.0" S="1.25"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s" name="distribution" x="@gammaShape.s:%s">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s" name="distr"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');

%                         fprintf(g, '\t\t\t\t\t\t\t\t<prior id="MutPrior.s:%s" name="distribution" x="@mutationRate.s:%s">\n', segments{s}, segments{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistribution.%s"  meanInRealSpace="true" M="0.0005" S="2" name="distr"/>\n', segments{s});
%                         fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    end
                end
                
                
                % constrain on the clades being monophyletic
                for i = 1:3
                    if sum(group(use_seqs)==i)>0
                        fprintf(g, '\t\t\t\t\t\t\t\t<distribution id="group%d" monophyletic="true"  spec="beast.math.distributions.MRCAPrior" tree="@%s.tree">\n', i, segments{1});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<taxonset id="group%d.set" spec="TaxonSet">\n',i);
                        for j = 1 : length(use_seqs)
                            if group(use_seqs(j))==i
                                fprintf(g, '\t\t\t\t\t\t\t\t\t\t<taxon idref="%s" spec="Taxon"/>\n', ids{use_seqs(j)});
                            end
                        end
                        fprintf(g, '\t\t\t\t\t\t\t\t\t</taxonset>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t\t</distribution>\n');
                    end
                end

            elseif contains(line, 'insert_operators')
                for s = 1 : length(segments)
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
%                                 fprintf(g, '\t\t\t\t<operator id="MutationScaler.s:%s.%d" spec="ScaleOperator" parameter="@mutationRate.s:%s.%d" scaleFactor="0.5" weight="1"/>\n', segments{s},i, segments{s},i);
                                fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s.%d" spec="ScaleOperator" parameter="@kappa.s:%s.%d" scaleFactor="0.5" weight="0.1"/>\n', segments{s},i, segments{s},i);
                                fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s.%d" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s.%d"/>\n', segments{s},i, segments{s},i);
                                fprintf(g, '\t\t\t\t<operator id="alpha_scaler.%s.%d" spec="ScaleOperator" parameter="@gammaShape.s:%s.%d" scaleFactor="0.75" weight="0.1"/>\n', segments{s},i, segments{s},i);

                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t<operator id="MutationScaler.s:%s" spec="ScaleOperator" parameter="@mutationRate.s:%s" scaleFactor="0.5" weight="1"/>\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s" spec="ScaleOperator" parameter="@kappa.s:%s" scaleFactor="0.5" weight="0.1"/>\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s"/>\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t<operator id="alpha_scaler.%s" spec="ScaleOperator" parameter="@gammaShape.s:%s" scaleFactor="0.75" weight="0.1"/>\n', segments{s}, segments{s});

                    end
                end
                
                
%                     fprintf(g, '\t\t\t\t<operator id="NetworkScaleRootOnly2" spec="NetworkScaleOperator" network="@network" scaleRootOnly="true" weight="3.0">\n');
%                     fprintf(g, '\t\t\t\t\t<downParameter idref="mutationRate.s:%s"/>\n', segments{1});
% 
%                     for i = 1 : length(segments)
%                          fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', segments{i});                 
%                     end
%                     fprintf(g, '\t\t\t\t</operator>\n');

             elseif contains(line, 'insert_mut_par')
                for s = 1 : length(segments)    
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
%                                 fprintf(g, '\t\t\t\t\t<downParameter idref="mutationRate.s:%s.%d"/>\n', segments{s},i);
                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t\t<downParameter idref="mutationRate.s:%s"/>\n', segments{s});
                    end
                end
             elseif contains(line, 'insert_muts_par')
                for s = 1 : length(segments)       
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
%                                 fprintf(g, '\t\t\t\t\t<down idref="mutationRate.s:%s.%d"/>\n', segments{s},i);
                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t\t<down idref="mutationRate.s:%s"/>\n', segments{s});
                    end
                end
             elseif contains(line, 'insert_param_log')
                for s = 1 : length(segments)                   
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
                                fprintf(g, '\t\t\t\t<log idref="treeLikelihood.%s.%d"/>\n', segments{s}, i);

                                fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s.%d"/>\n', segments{s}, i);      
                                fprintf(g, '\t\t\t\t<log idref="kappa.s:%s.%d"/>\n', segments{s},i);
                                fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s.%d"/>\n', segments{s},i);
                                fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s.%d"/>\n', segments{s},i);

                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t<log idref="treeLikelihood.%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="kappa.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s"/>\n', segments{s});

                    end

                end
            elseif contains(line, 'insert_tree_likelihood')
                for s = 1 : length(segments)
                    if s==1 
                        for i = 1:3
                            if sum(group(use_seqs)==i)>0
                                if i==1
                                    mutrate= length(fasta3(1).Sequence)/sum(~const1);
                                    range=[1, sum(~const1)];
                                elseif i==3
                                    mutrate=1;
                                    range=[sum(~const1)+sum(~const2)+1, length(fasta_comb(1).Sequence)-10];
                                else
                                    er
                                end
                                
                                fprintf(g, '\t\t\t\t\t<distribution id="treeLikelihood.%s.%d" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',segments{s},i,segments{s});
                                fprintf(g, '\t\t\t\t\t\t<data id="%s.%d" spec="FilteredAlignment" filter="%d-%d" data="@%s"/>\n', segments{s}, i, range(1), range(2), segments{s});
                                fprintf(g, '\t\t\t\t\t\t<siteModel id="SiteModel.s:%s.%d" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s.%d">\n',segments{s},i,segments{s},i);
%                                 fprintf(g, '\t\t\t\t\t\t<parameter spec="parameter.RealParameter" estimate="false" name="mutationRate">%f</parameter>\n', mutrate);

                                fprintf(g, '\t\t\t\t\t\t\t<substModel id="hky.s:%s.%d" spec="HKY" kappa="@kappa.s:%s.%d">\n',segments{s},i,segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s.%d" spec="Frequencies" frequencies="@freqParameter.s:%s.%d"/>\n',segments{s},i,segments{s},i);
                                fprintf(g, '\t\t\t\t\t\t\t</substModel>\n');
                                fprintf(g, '\t\t\t\t\t\t</siteModel>\n');
                                fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s.%d" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@mutationRate.s:%s.%d"/>\n',segments{s},i,segments{s},i);
                                fprintf(g, '\t\t\t\t\t</distribution>\n');
                            end
                        end
                    else
                        fprintf(g, '\t\t\t\t\t<distribution id="treeLikelihood.%s" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true" data="@%s">\n',segments{s},segments{s},segments{s});
                        fprintf(g, '\t\t\t\t\t\t<siteModel id="SiteModel.s:%s" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s">\n',segments{s},segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.s:%s">\n',segments{s},segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s" spec="Frequencies" frequencies="@freqParameter.s:%s"/>\n',segments{s},segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t</substModel>\n');
                        fprintf(g, '\t\t\t\t\t\t</siteModel>\n');
                        fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@mutationRate.s:%s"/>\n',segments{s},segments{s});
                        fprintf(g, '\t\t\t\t\t</distribution>\n');
                    end
                end
            elseif contains(line, 'insert_weights')
                for i = 1 : length(segments)
                    if startsWith(segs_files(i).name,'LN')
                        fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(i).name]);
                    else
                        fasta = fastaread(['data/Aln_Trees/FlexSon_Ss046/' segs_files(i).name]);
                    end

                    fprintf(g, '%d ', length(fasta(1).Sequence));
                end
            elseif contains(line, 'insert_seg_tree_state')
                for i = 1 : length(segments)
                    fprintf(g, '\t\t\t\t<tree id="%s.tree" spec="beast.evolution.tree.Tree" name="stateNode">\n', segments{i});

                    fprintf(g, '\t\t\t\t\t<taxonset id="TaxonSet.%s" spec="TaxonSet">\n', segments{i});
                    fprintf(g, '\t\t\t\t\t\t<alignment idref="%s"/>\n', segments{i});
                    fprintf(g, '\t\t\t\t\t</taxonset>\n');
                    fprintf(g, '\t\t\t\t\t<trait idref="traitSet"/>\n');

                    fprintf(g, '\t\t\t\t</tree>\n');
                end

            elseif contains(line, 'insert_init_trees')
                for i = 1 : length(segments)
                     fprintf(g, '\t\t\t\t\t<tree idref="%s.tree"/>\n', segments{i});                 
                end

            elseif contains(line, 'insert_seg_tree')
                for i = 1 : length(segments)
                     fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', segments{i});                 
                end
            elseif contains(line, 'insert_seg_logger')
                for i = 1 : length(segments)
                    fprintf(g, '\t\t\t<logger spec="Logger" logEvery="50000" mode="tree" fileName="$(filebase).%s.trees">\n', segments{i});           
                    fprintf(g, '\t\t\t\t<log idref="%s.tree"/>\n', segments{i});           
                    fprintf(g, '\t\t\t</logger>\n');
                end

            elseif contains(line, 'insert_stats_log')
                for i = 1 : length(segments)
                    fprintf(g, '\t\t\t\t<log spec="TreeStatLogger" tree="@%s.tree"/>\n', segments{i});      
                end
            elseif contains(line, 'insert_seg_tree')
                for i = 1 : length(segments)
                     fprintf(g, '\t\t\t\t\t<segmentTree idref="%s.tree"/>\n', segments{i});                 
                end

            else
                fprintf(g, line);
            end
        end
        fclose(f);fclose(g);

    end
    
    
    
    %% make 3 replicates of the tree inference xmls    
    for r = 0 : 2
        for seg = 1 :length(segments)
            % build the inference xml
            f_inf = fopen('individualtrees_template.xml');

            % open the inference xml
            g = fopen(sprintf('xmls/SonFlex_%s_rep%d.xml', segments{seg}, r), 'w');

            % keep track of the segment count for the nexus file name
            segmentcount = 1; segmentcount2=1;

            while ~feof(f_inf)
                line = fgets(f_inf);
                if ~isempty(strfind(line, 'insert_taxa'))
                elseif contains(line, 'insert_sequence')
                    fprintf(g, '\t<data id="alignment">\n');
                    if startsWith(segs_files(seg).name,'LN')
                        fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(seg).name]);
                    else
                        fasta = fastaread(['data/Aln_Trees/FlexSon_Ss046/' segs_files(seg).name]);       
                    end

                    seq_length(seg) = length(fasta(1).Sequence);

                    for j = 1 : length(use_seqs)
                        ind = has_plasmid(use_seqs(j),seg);                        
                        if ind~=-1 && group(use_seqs(j))==1
                            fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                 segments{seg}, fasta(ind).Header,...
                                 fasta(ind).Header, fasta(ind).Sequence);
                        end
                    end
                    fprintf(g, '\t</data>\n');
                elseif ~isempty(strfind(line, 'insert_sampling_times'))
                    indices = use_seqs(has_plasmid(use_seqs,seg)>-1 & group(use_seqs)'==1);

                    for j = 1 : length(indices)

                        if j==length(indices)
                            fprintf(g,'\t\t\t\t%s=%s-01\n', ids{indices(j)}, times{indices(j)});
                        else
                            fprintf(g,'\t\t\t\t%s=%s-01,\n', ids{indices(j)}, times{indices(j)});
                        end
                    end 
                elseif contains(line, 'id="clockRate')
                    fprintf(g, strrep(line, '0.0001', num2str(clockrates(seg))));
                elseif contains(line, 'insert_weights')
                    fprintf(g, strrep(line, 'insert_weights', num2str(slen)));
                elseif ~isempty(strfind(line, 'initial_Ne'))
                    fprintf(g, strrep(line, 'initial_Ne', num2str(exprnd(1))));
                elseif contains(line, 'insert_nexus_name')
                    line = strrep(line, 'insert_nexus_name', sprintf('sim_%d.%s.alignment.nexus', i,segment_names{seg}) );
                    fprintf(g, line);
%                 elseif contains(line, 'ClockRateOperator')
%                     line = fgets(f);
%                     line = fgets(f);
%                 elseif contains(line, '<up idref="clockRate"/>')
                    
                else
                    fprintf(g, '%s', line);
                end
            end

            fclose(f_inf); fclose(g);
        end
    end

end


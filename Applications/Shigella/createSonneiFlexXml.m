rng(187989)

clear

weights_name = {'equal','weighted'};


clockrates = [ 9.1586e-05 10000000 3.3748e-04];
for w = 1

    for cutoff = 50

        %% Shigella sonnei
        meta1 = importdata('data/Sonnei_MonthYear_07082022.txt');
        meta2 = importdata('data/ShigFlex_01092022_YearMonth.txt');

        flex = fastaread('data/Flex/Sflex_NC_004337_21022024_NoSero6_cleanGubbinsV241.filtered_polymorphic_sites.fasta');
        sonnei = fastaread('data/Sonnei/Sson_NC_007384_18122023_cleanGubbinsV241.filtered_polymorphic_sites.snpsites.fasta');

        % Create a combine fasta file
        c=1;    
        clear fasta_comb;    
        for i = 1 : length(flex)
            fasta_comb(c).Header = flex(i).Header;
            fasta_comb(c).Sequence = [flex(i).Sequence, repmat('N', 1,length(sonnei(1).Sequence))];
            c=c+1;
        end
        for i = 1 : length(sonnei)
            fasta_comb(c).Header = sonnei(i).Header;
            fasta_comb(c).Sequence = [repmat('N', 1,length(flex(1).Sequence)) sonnei(i).Sequence];
            c=c+1;
        end

        delete('data/chromosome_concat.aln');
        fastawrite('data/chromosome_concat.aln',fasta_comb);

        NC_007384 = struct2table(fastaread('data/chromosome_concat.aln'));
        LN624486 = struct2table(fastaread('data/pksr100/SonFlex_LN624486_over80.filtered95.aln'));

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
                    error('should not occur')
                end
            end
            if ismember(ids{i}, LN624486.Header)
                has_plasmid(i,2) = find(ismember(LN624486.Header,ids{i}));
            end
    %         if ismember(ids{i}, NC_007385.Header)
    %             has_plasmid(i,3) = find(ismember(NC_007385.Header,ids{i}));
    %         end

        end
    
        if w==2
            % only use sequences if they have the psk100 plasmid
            dont_use_seqs = ~(has_plasmid(:,2)>-1 & group'~=2);
        else
            dont_use_seqs = ~(has_plasmid(:,2)>-2 & group'~=2);
        end

        plasmid_group = group;
        plasmid_group(dont_use_seqs) =-1;

        for i = 1 : max(plasmid_group)
            weights(plasmid_group==i) = 1/sum(plasmid_group==i);
        end
        weights(dont_use_seqs) = 0;

        indices = [];
        for i = 1 : 400
            use_weights = weights;
            use_weights(indices) = 0;
            indices = [indices, randsample(length(weights), 1, true, use_weights)];
        end
        use_seqs = sort(indices);

        segs_files(1).name = 'chromosome_concat.aln';
        segs_files(2).name = 'SonFlex_LN624486_over80.filtered95.aln'
    %     segs_files(3).name = 'SonFlex_NC_007385_Over70_23082022.aln';

        segments = {'NC_007384','LN624486'};

    %     dsa

        %%
        for r = 0:2
            f = fopen('template.xml');
            g = fopen(['xmls/SonFlex_' weights_name{w} '_rep' num2str(r) '.xml'], 'w');

            est_tip_time = cell(0,0);
            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_alignment')
                    for i = 1 : length(segments)
                        fprintf(g, '\t<data id="%s">\n',segments{i});
                        if startsWith(segs_files(i).name,'SonFlex')
                            fasta = fastaread(['data/pksr100/' segs_files(i).name]);
                        else
                            fasta = fastaread(['data/' segs_files(i).name]);
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
                                length_vals = length(sonnei(1).Sequence);
                                start_vals = length(flex(1).Sequence)+1;

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
                    fprintf(g, '\t\t<run id="mcmc" spec="coupledMCMC.CoupledMCMC" chainLength="10000000" storeEvery="1000000">\n');            
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
                elseif contains(line, 'id="clockRate.c" name="stateNode">0.0005</parameter>')
                    fprintf(g, strrep(line, '0.0005', '1.0'));
                elseif contains(line, 'insert_parameters')

                    for s = 1 : length(segments)
                        if s==1 
    %                         fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" name="stateNode" lower="2.5e-04" upper="4e-04">3.5e-04</parameter>\n',segments{s});

                            for i = 1:3
                                if sum(group(use_seqs)==i)>0
                                    if i==3
                                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="mutationRate.s:%s.%d" name="stateNode">%f</parameter>\n',segments{s}, i, clockrates(i));
                                    else
                                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="mutationRate.s:%s.%d" name="stateNode">%f</parameter>\n',segments{s}, i, clockrates(i));
                                    end
                                    fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="kappa.s:%s.%d" lower="0.0" name="stateNode">%f</parameter>\n',segments{s},i, lognrnd(0,0.5,1));
                                    fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="gammaShape.s:%s.%d" name="stateNode">%f</parameter>\n',segments{s},i, lognrnd(0,0.5,1));
                                    fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="freqParameter.s:%s.%d" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',segments{s},i);

                                end
                            end
                        else
                            fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="mutationRate.s:%s" name="stateNode">1</parameter>\n',segments{s});
                            fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="kappa.s:%s" lower="0.0" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                            fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="gammaShape.s:%s" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                            fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="freqParameter.s:%s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',segments{s});
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
                            fprintf(g, '\t\t\t\t\t\t\t\t<distribution id="group%d" monophyletic="true"  spec="beast.base.evolution.tree.MRCAPrior" tree="@%s.tree">\n', i, segments{1});
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
                    fprintf(g, '\t\t\t<operator id="AVMNOperator" spec="kernel.AdaptableVarianceMultivariateNormalOperator" allowNonsense="true" beta="0.05" burnin="400" initial="800" weight="0.01">\n');
                    fprintf(g, '\t\t\t\t<transformations id="AVMNSumTradnsform.h1n1pdm_NA" spec="operator.kernel.Transform$LogConstrainedSumTransform">\n');
                    for s = 1 : length(segments)
                        if s==1 
                            for i = 1:3
                                if sum(group(use_seqs)==i)>0                                
                                    fprintf(g, '\t\t\t\t\t<f idref="freqParameter.s:%s.%d"/>\n', segments{s},i);
                                end
                            end
                        else
                            fprintf(g, '\t\t\t\t\t<f idref="freqParameter.s:%s"/>\n', segments{s});
                        end
                    end
                    fprintf(g, '\t\t\t\t</transformations>\n');
                    fprintf(g, '\t\t\t\t<transformations id="AVMNLogTransform.h1n1pdm_NA" spec="operator.kernel.Transform$LogTransform">\n');
                    for s = 1 : length(segments)
                        % fprintf(g, '\t\t\t\t\t<f idref="clockRate.c"/>\n');
                        if s==1 
                            for i = 1:3
                                if sum(group(use_seqs)==i)>0
                                    fprintf(g, '\t\t\t\t\t<f idref="kappa.s:%s.%d"/>\n', segments{s},i);
                                    fprintf(g, '\t\t\t\t\t<f idref="gammaShape.s:%s.%d"/>\n', segments{s},i);
                                end
                            end
                        else
                            fprintf(g, '\t\t\t\t\t<f idref="mutationRate.s:%s"/>\n', segments{s});
                            fprintf(g, '\t\t\t\t\t<f idref="kappa.s:%s"/>\n', segments{s});
                            fprintf(g, '\t\t\t\t\t<f idref="gammaShape.s:%s"/>\n', segments{s});
                        end
                    end
                    fprintf(g, '\t\t\t\t</transformations>\n');
                    fprintf(g, '\t\t\t\t<kernelDistribution id="KernelDistribution$Bactrian.7" spec="operator.kernel.KernelDistribution$Bactrian"/>\n');                
                    fprintf(g, '\t\t\t</operator>\n');

                    % fprintf(g, '\t\t\t\t<operator id="ClockRateScaler.s" spec="AdaptableOperatorSampler" parameter="@clockRate.c" weight="0.1" operator="@AVMNOperator">\n');
                    % fprintf(g, '\t\t\t\t\t<operator id="ClockRateScalerX.s" spec="kernel.BactrianScaleOperator" parameter="@clockRate.c" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n');
                    % fprintf(g, '\t\t\t\t</operator>\n');


                    for s = 1 : length(segments)
                        if s==1 
                            for i = 1:3
                                if sum(group(use_seqs)==i)>0            
%                                     fprintf(g, '\t\t\t\t<operator id="MutationScaler.s:%s.%d" spec="AdaptableOperatorSampler" parameter="@mutationRate.s:%s.%d" weight="1" operator="@AVMNOperator">\n', segments{s},i, segments{s},i);
%                                     fprintf(g, '\t\t\t\t\t<operator id="MutationScalerX.s:%s.%d" spec="kernel.BactrianScaleOperator" parameter="@mutationRate.s:%s.%d" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},i,segments{s},i);
%                                     fprintf(g, '\t\t\t\t</operator>\n');

                                    fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s.%d" spec="AdaptableOperatorSampler" parameter="@kappa.s:%s.%d" weight="0.01" operator="@AVMNOperator">\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t\t<operator id="KappaScalerX.s:%s.%d" spec="kernel.BactrianScaleOperator" parameter="@kappa.s:%s.%d" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t</operator>\n');

                                    fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s.%d" spec="AdaptableOperatorSampler" weight="0.01" parameter="@freqParameter.s:%s.%d" operator="@AVMNOperator">\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t\t<operator id="FrequenciesExchangerX.s:%s.%d" spec="operator.kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s.%d"/>\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t</operator>\n');


                                    fprintf(g, '\t\t\t\t<operator id="alpha_scaler.%s.%d" spec="AdaptableOperatorSampler" parameter="@gammaShape.s:%s.%d" weight="0.01" operator="@AVMNOperator">\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t\t<operator id="alpha_scalerX.s:%s.%d" spec="kernel.BactrianScaleOperator" parameter="@gammaShape.s:%s.%d" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},i, segments{s},i);
                                    fprintf(g, '\t\t\t\t</operator>\n');


                                end
                            end
                        else

                        fprintf(g, '\t\t\t\t<operator id="MutationRateScaler.s:%s" spec="AdaptableOperatorSampler" parameter="@mutationRate.s:%s" weight="0.01" operator="@AVMNOperator">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t<operator id="MutationRateScalerX.s:%s" spec="kernel.BactrianScaleOperator" parameter="@mutationRate.s:%s" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},segments{s});
                        fprintf(g, '\t\t\t\t</operator>\n');


                        fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s" spec="AdaptableOperatorSampler" parameter="@kappa.s:%s" weight="0.01" operator="@AVMNOperator">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t<operator id="KappaScalerX.s:%s" spec="kernel.BactrianScaleOperator" parameter="@kappa.s:%s" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},segments{s});
                        fprintf(g, '\t\t\t\t</operator>\n');

                        fprintf(g, '\t\t\t\t<operator id="FrequenciesExchanger.s:%s" spec="AdaptableOperatorSampler" weight="0.01" parameter="@freqParameter.s:%s" operator="@AVMNOperator">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t<operator id="FrequenciesExchangerX.s:%s" spec="operator.kernel.BactrianDeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:%s"/>\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t</operator>\n');

                        fprintf(g, '\t\t\t\t<operator id="alpha_scaler.%s" spec="AdaptableOperatorSampler" parameter="@gammaShape.s:%s" weight="0.01" operator="@AVMNOperator">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t<operator id="alpha_scalerX.s:%s" spec="kernel.BactrianScaleOperator" parameter="@gammaShape.s:%s" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},segments{s});
                        fprintf(g, '\t\t\t\t</operator>\n');


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
                                    fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s.%d"/>\n', segments{s},i);
                                end
                            end
                        else
                            fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s"/>\n', segments{s});
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
                    length_for_weights = zeros(0,0);
                    for s = 1 : length(segments)
                        if s==1 
                            for i = 1:3
                                if sum(group(use_seqs)==i)>0
                                    if i==1
                                        range=[1, length(flex(1).Sequence)];
                                    elseif i==3
                                        range=[length(flex(1).Sequence)+1, length(fasta_comb(1).Sequence)];
                                    end

                                    length_for_weights = [length_for_weights range(2)-range(1)];

                                    fprintf(g, '\t\t\t\t\t<distribution id="treeLikelihood.%s.%d" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true">\n',segments{s},i,segments{s});
                                    fprintf(g, '\t\t\t\t\t\t<data id="%s.%d" spec="FilteredAlignment" filter="%d-%d" data="@%s"/>\n', segments{s}, i, range(1), range(2), segments{s});
                                    fprintf(g, '\t\t\t\t\t\t<siteModel id="SiteModel.s:%s.%d" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s.%d" mutationRate="@mutationRate.s:%s.%d">\n',segments{s},i,segments{s},i,segments{s},i);
    %                                 fprintf(g, '\t\t\t\t\t\t<parameter spec="parameter.RealParameter" estimate="false" name="mutationRate">%f</parameter>\n', mutrate);

                                    fprintf(g, '\t\t\t\t\t\t\t<substModel id="hky.s:%s.%d" spec="HKY" kappa="@kappa.s:%s.%d">\n',segments{s},i,segments{s},i);
                                    fprintf(g, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s.%d" spec="Frequencies" frequencies="@freqParameter.s:%s.%d"/>\n',segments{s},i,segments{s},i);
                                    fprintf(g, '\t\t\t\t\t\t\t</substModel>\n');
                                    fprintf(g, '\t\t\t\t\t\t</siteModel>\n');
                                    fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s.%d" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',segments{s},i);
                                    fprintf(g, '\t\t\t\t\t</distribution>\n');
                                end
                            end
                        else
                            fprintf(g, '\t\t\t\t\t<distribution id="treeLikelihood.%s" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true" data="@%s">\n',segments{s},segments{s},segments{s});
                            fprintf(g, '\t\t\t\t\t\t<siteModel id="SiteModel.s:%s" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s" mutationRate="@mutationRate.s:%s">\n',segments{s},segments{s},segments{s});
                            fprintf(g, '\t\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.s:%s">\n',segments{s},segments{s});
                            fprintf(g, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s" spec="Frequencies" frequencies="@freqParameter.s:%s"/>\n',segments{s},segments{s});
                            fprintf(g, '\t\t\t\t\t\t\t</substModel>\n');
                            fprintf(g, '\t\t\t\t\t\t</siteModel>\n');
                            fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c"/>\n',segments{s});
                            fprintf(g, '\t\t\t\t\t</distribution>\n');
                        end
                    end
                elseif contains(line, 'insert_weights')
                    fprintf(g, '%d %d ', length_for_weights(1), length_for_weights(2));
                    for i = 1 : length(segments)
                        if startsWith(segs_files(i).name,'LN')
                            fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(i).name]);
                            fprintf(g, '%d ', length(fasta(1).Sequence));
                        end
                    end
                    fprintf(g, '\n');

                elseif contains(line, 'insert_seg_tree_state')
                    for i = 1 : length(segments)
                        fprintf(g, '\t\t\t\t<tree id="%s.tree" spec="beast.base.evolution.tree.Tree" name="stateNode">\n', segments{i});

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
                elseif contains(line, '<downParameter idref="clockRate.c"/>')
                    % skip this line
                elseif contains(line, '<operator id="FixMeanMutationRatesOperator"')
                    fgets(f);fgets(f);fgets(f);fgets(f);fgets(f);
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



%         %% make 3 replicates of the tree inference xmls    
%         for r = 0 : 2
%             for seg = 1 :length(segments)
%                 % build the inference xml
%                 f_inf = fopen('individualtrees_template.xml');
% 
%                 % open the inference xml
%                 g = fopen(sprintf('xmls/SonFlex_%s_%s_rep%d.xml',weights_name{w}, segments{seg}, r), 'w');
% 
%                 % keep track of the segment count for the nexus file name
%                 segmentcount = 1; segmentcount2=1;
% 
%                 while ~feof(f_inf)
%                     line = fgets(f_inf);
%                     if ~isempty(strfind(line, 'insert_taxa'))
%                     elseif contains(line, 'insert_sequence')
%                         if startsWith(segs_files(seg).name,'LN')
%                             fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(seg).name]);
%                         else
%                             fasta = fastaread(['data/Aln_Trees/FlexSon_Ss046/' segs_files(seg).name]);       
%                         end
% 
%                         seq_length(seg) = length(fasta(1).Sequence);
% 
%                         for j = 1 : length(use_seqs)
%                             ind = has_plasmid(use_seqs(j),seg);                        
%                             if ind~=-1 && group(use_seqs(j))==1
%                                 fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
%                                      segments{seg}, fasta(ind).Header,...
%                                      fasta(ind).Header, fasta(ind).Sequence);
%                             end
%                         end
%                     elseif ~isempty(strfind(line, 'insert_sampling_times'))
%                         indices = use_seqs(has_plasmid(use_seqs,seg)>-1);
%                         stimes = 'rem';
% 
%                         for j = 1 : length(indices)
%                             stimes = [stimes ',' ids{indices(j)} '=' times{indices(j)} '-01'];
%                         end 
%                         fprintf(g, strrep(line, 'insert_sampling_times', strrep(stimes, 'rem,','')));
%                     elseif contains(line, 'id="clockRate')
%                         fprintf(g, strrep(line, '0.0001', num2str(clockrates(seg))));
%                     elseif contains(line, 'insert_weights')
%                         fprintf(g, strrep(line, 'insert_weights', num2str(slen)));
%                     elseif ~isempty(strfind(line, 'initial_Ne'))
%                         fprintf(g, strrep(line, 'initial_Ne', num2str(exprnd(1))));
%                     elseif contains(line, 'insert_nexus_name')
%                         line = strrep(line, 'insert_nexus_name', sprintf('sim_%d.%s.alignment.nexus', i,segment_names{seg}) );
%                         fprintf(g, line);
%     %                 elseif contains(line, 'ClockRateOperator')
%     %                     line = fgets(f);
%     %                     line = fgets(f);
%     %                 elseif contains(line, '<up idref="clockRate"/>')
% 
%                     else
%                         fprintf(g, '%s', line);
%                     end
%                 end
% 
%                 fclose(f_inf); fclose(g);
%             end
%         end

    end
end


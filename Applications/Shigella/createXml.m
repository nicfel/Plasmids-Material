rng(187989)

clear

weights_name = {'equal','weighted'};

for w = 1 : 1
    for dataset_nr = 1 : 2

        if dataset_nr==1
            %% Shigella sonnei
            meta = importdata('data/Sonnei_MonthYear_07082022.txt');


            NC_007384 = struct2table(fastaread('data/Aln_Trees/Ss046_Sonnei_ref/Sson_NC_007384_29042022clean_gubbinsv241.filtered_polymorphic_sites.fasta'));
            NC_007385 = struct2table(fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_007385_plasmid_SonneiOnly_Over70_0408202.aln'));
            NC_009345 = struct2table(fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009345_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln'));
            NC_009346 = struct2table(fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009346_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln'));
            NC_009347 = struct2table(fastaread('data/Aln_Trees/Ss046_Sonnei_ref/NC_009347_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln'));
            ids = meta.textdata(2:end,1);
            times = strrep(meta.textdata(2:end,2), '_','-');


            has_plasmid = -1*ones(length(ids), 6);
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

            weights = (has_plasmid(:,2)>0) + (has_plasmid(:,3)>0) + (has_plasmid(:,4)>0) + (has_plasmid(:,5)>0)+ (has_plasmid(:,6)>0);
            
            if w==1
                weights = has_plasmid(:,1)>0;                
            end
            
            indices = [];
            for i = 1 : 200
                use_weights = weights;
                use_weights(indices) = 0;
                indices = [indices, randsample(length(weights), 1, true, use_weights)];
            end
            use_seqs = sort(indices);


            segs_files(1).name = 'Sson_NC_007384_29042022clean_gubbinsv241.filtered_polymorphic_sites.fasta';
            segs_files(2).name = 'NC_007385_plasmid_SonneiOnly_Over70_0408202.aln';
            segs_files(3).name = 'NC_009345_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln';
            segs_files(4).name = 'NC_009346_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln';
            segs_files(5).name = 'NC_009347_plasmid_SonneiOnly_Over50_26082022_full_clean30092022.aln';
            segments = {'NC_007384','NC_007385','NC_009345', 'NC_009346', 'NC_009347'};

            f = fopen('sonnei_plasmid_info.tsv','w');
            fprintf(f, 'id\ttimes\t%s\t%s\t%s\t%s\t%s\n',segments{1},segments{2},segments{3},segments{4},segments{5});
            for i = 1 : size(has_plasmid,1)
                fprintf(f, '%s\t%s-01\t%d\t%d\t%d\t%d\t%d\n',ids{i},times{i},...
                    has_plasmid(i,1)~=-1,has_plasmid(i,2)~=-1,has_plasmid(i,3)~=-1,has_plasmid(i,4)~=-1,has_plasmid(i,5)~=-1);
            end
            fclose(f);


        elseif dataset_nr==2
            %% Shigella felxneri        
            meta = importdata('data/ShigFlex_01092022_YearMonth.txt');                        
            NC_004337 = struct2table(fastaread('data/Aln_Trees/Flexneri_ref/Flex1_NC_007384.aln'));
            NC_004851 = struct2table(fastaread('data/Aln_Trees/Flexneri_ref/NC_004851_plasmid_FlexneriOnly_Over70_24082022.aln'));
            NC_004851 = struct2table(fastaread('data/Aln_Trees/Flexneri_ref/NC_004851_plasmid_FlexneriOnly_Over70_24082022.aln'));


            ids = meta.textdata(2:end,1);
            times = strrep(meta.textdata(2:end,2), '_','-');

            has_plasmid = -1*ones(length(ids), 6);
            for i = 1 : length(ids)
                if ismember(ids{i}, NC_004337.Header)
                    has_plasmid(i,1) = find(ismember(NC_004337.Header,ids{i}));
                end
                if ismember(ids{i}, NC_004851.Header)
                    has_plasmid(i,2) = find(ismember(NC_004851.Header,ids{i}));
                end
            end

            use_seqs = find(has_plasmid(:,1)>0);
            use_seqs = sort(randsample(use_seqs, min(length(use_seqs),200)));

            segs_files(1).name = 'Flex1_NC_007384.aln';
            segs_files(2).name = 'NC_004851_plasmid_FlexneriOnly_Over70_24082022.aln';
            segments = {'NC_004337','NC_004851'};


        end

        %%
        for r = 0:2
            f = fopen('template.xml');
            if dataset_nr==1                
                g = fopen(['xmls/Sonnei_' weights_name{w} '_rep' num2str(r) '.xml'], 'w');
            else
                g = fopen(['xmls/Flexneri_' weights_name{w} '_rep' num2str(r) '.xml'], 'w');
            end

            est_tip_time = cell(0,0);
            while ~feof(f)
                line = fgets(f);
                if contains(line, 'insert_alignment')
                    for i = 1 : length(segments)
                        fprintf(g, '\t<data id="%s">\n',segments{i});
                        if startsWith(segs_files(i).name,'LN')
                            fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(i).name]);
                        elseif dataset_nr==1
                            fasta = fastaread(['data/Aln_Trees/Ss046_Sonnei_ref/' segs_files(i).name]);
                        else
                            fasta = fastaread(['data/Aln_Trees/Flexneri_ref/' segs_files(i).name]);
                        end

                        seq_length(i) = length(fasta(1).Sequence);

                        for j = 1 : length(use_seqs)
                            ind = has_plasmid(use_seqs(j),i);
                            if ind==-1
        %                         fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
        %                              segments{i}, ids{use_seqs(j)},...
        %                              ids{use_seqs(j)}, repmat('N',1,length(fasta(1).Sequence)));
                            else
                                fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                     segments{i}, fasta(ind).Header,...
                                     fasta(ind).Header, fasta(ind).Sequence);
                            end
                        end
                        fprintf(g, '\t</data>\n');

                    end
                elseif contains(line, 'insert_run_header')
                 fprintf(g, '\t\t<run id="mcmc" spec="coupledMCMC.CoupledMCMC" chainLength="10000000" deltaTemperature="0.00001" heatLikelihoodOnly="true" storeEvery="1000000" chains="2" resampleEvery="10000">\n');            
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
                elseif contains(line, 'id="plasmidTransferRate"')
                    fprintf(g, strrep(line, 'id="plasmidTransferRate"', sprintf('id="plasmidTransferRate" dimension="%d"', length(segments)-1)));
                elseif contains(line, 'insert_nr_segments')
                    fprintf(g, strrep(line, 'insert_nr_segments', num2str(length(segments))));
                elseif contains(line, 'insert_parameters')

                    for s = 1 : length(segments)
                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="kappa.s:%s" lower="0.0" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="mutationRate.s:%s" name="stateNode">0.0005</parameter>\n',segments{s});
                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="gammaShape.s:%s" name="stateNode">%f</parameter>\n',segments{s}, lognrnd(0,0.5,1));
                        fprintf(g, '\t\t\t\t<parameter spec="parameter.RealParameter" id="freqParameter.s:%s" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n',segments{s});
                    end
                elseif contains(line, 'insert_corename')
                    fprintf(g, strrep(line, 'insert_corename', segments{1}));
                elseif contains(line, 'insert_priors')
                    % insert sampling time priors
                    for s = 1 : length(segments)
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="KappaPrior.s:%s" name="distribution" x="@kappa.s:%s">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.%s.1" name="distr" M="1.0" S="1.25"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="GammaPrior.s:%s" name="distribution" x="@gammaShape.s:%s">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<Exponential id="ExponentialDistribution.%s" name="distr"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                        fprintf(g, '\t\t\t\t\t\t\t\t<prior id="MutPrior.s:%s" name="distribution" x="@mutationRate.s:%s">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t\t<LogNormal id="LogNormalDistribution.%s"  meanInRealSpace="true" M="0.0005" S="2" name="distr"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t\t\t\t</prior>\n');
                    end
                elseif contains(line, 'insert_operators')
                    fprintf(g, '\t\t\t<operator id="AVMNOperator" spec="kernel.AdaptableVarianceMultivariateNormalOperator" allowNonsense="true" beta="0.05" burnin="400" initial="800" weight="0.01">\n');
                    fprintf(g, '\t\t\t\t<transformations id="AVMNSumTradnsform.h1n1pdm_NA" spec="operator.kernel.Transform$LogConstrainedSumTransform">\n');
                    for s = 1 : length(segments)
                        fprintf(g, '\t\t\t\t\t<f idref="freqParameter.s:%s"/>\n', segments{s});
                    end
                    fprintf(g, '\t\t\t\t</transformations>\n');
                    fprintf(g, '\t\t\t\t<transformations id="AVMNLogTransform.h1n1pdm_NA" spec="operator.kernel.Transform$LogTransform">\n');
                    for s = 1 : length(segments)
                        fprintf(g, '\t\t\t\t\t<f idref="mutationRate.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t<f idref="kappa.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t\t<f idref="gammaShape.s:%s"/>\n', segments{s});
                    end
                    fprintf(g, '\t\t\t\t</transformations>\n');
                    fprintf(g, '\t\t\t\t<kernelDistribution id="KernelDistribution$Bactrian.7" spec="operator.kernel.KernelDistribution$Bactrian"/>\n');                
                    fprintf(g, '\t\t\t</operator>\n');                
                    for s = 1 : length(segments)
                        fprintf(g, '\t\t\t\t<operator id="MutationScaler.s:%s" spec="AdaptableOperatorSampler" parameter="@mutationRate.s:%s" weight="0.1" operator="@AVMNOperator">\n', segments{s}, segments{s});
                        fprintf(g, '\t\t\t\t\t<operator id="MutationScalerX.s:%s" spec="kernel.BactrianScaleOperator" parameter="@mutationRate.s:%s" scaleFactor="0.1" upper="10.0" weight="0.1"/>\n', segments{s},segments{s});
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

                 elseif contains(line, 'insert_mut_par')
                     for s = 1 : length(segments)                   
                         fprintf(g, '\t\t\t\t\t<downParameter idref="mutationRate.s:%s"/>\n', segments{s});
                     end
                 elseif contains(line, 'insert_muts_par')
    %                 for s = 1 : length(segments)                   
    %                     fprintf(g, '\t\t\t\t\t<down idref="mutationRate.s:%s"/>\n', segments{s});
    %                 end
                 elseif contains(line, 'insert_param_log')
                    for s = 1 : length(segments)                   
                        fprintf(g, '\t\t\t\t<log idref="treeLikelihood.%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="kappa.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="mutationRate.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="gammaShape.s:%s"/>\n', segments{s});
                        fprintf(g, '\t\t\t\t<log idref="freqParameter.s:%s"/>\n', segments{s});
                    end
                elseif contains(line, 'insert_tree_likelihood')
                    for i = 1 : length(segments)
                        fprintf(g, '\t\t\t\t\t<distribution id="treeLikelihood.%s" spec="ThreadedTreeLikelihood" tree="@%s.tree" useAmbiguities="true" data="@%s">\n',segments{i},segments{i},segments{i});

                        fprintf(g, '\t\t\t\t\t\t<siteModel id="SiteModel.s:%s" spec="SiteModel"  gammaCategoryCount="4" shape="@gammaShape.s:%s">\n',segments{i},segments{i});
                        fprintf(g, '\t\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.s:%s">\n',segments{i},segments{i});
                        fprintf(g, '\t\t\t\t\t\t\t\t<frequencies id="estimatedFreqs.s:%s" spec="Frequencies" frequencies="@freqParameter.s:%s"/>\n',segments{i},segments{i});
                        fprintf(g, '\t\t\t\t\t\t\t</substModel>\n');
                        fprintf(g, '\t\t\t\t\t\t</siteModel>\n');
                        fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@mutationRate.s:%s"/>\n',segments{i},segments{i});

    %                         fprintf(g, '\t\t\t\t\t\t<branchRateModel id="StrictClock.%s" spec="plasmids.ratemodel.EarlyLateRate" finalTime="100" clock.rate="@mutationRate.s:%s">\n',segments{i},segments{i});
    %                         fprintf(g, '\t\t\t\t\t\t\t<parameter name="lateRate">0.001</parameter>\n');
    %                         fprintf(g, '\t\t\t\t\t\t</branchRateModel>\n');
                        fprintf(g, '\t\t\t\t\t</distribution>\n');
                    end
                elseif contains(line, 'insert_weights')
                    for i = 1 : length(segments)
                        if startsWith(segs_files(i).name,'LN')
                            fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(i).name]);
                        elseif dataset_nr==1
                            fasta = fastaread(['data/Aln_Trees/Ss046_Sonnei_ref/' segs_files(i).name]);
                        else
                            fasta = fastaread(['data/Aln_Trees/Flexneri_ref/' segs_files(i).name]);
                        end

                        fprintf(g, '%d ', length(fasta(1).Sequence));
                    end
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
                if dataset_nr==1
                    g = fopen(sprintf('xmls/Sonnei_%s_%s_rep%d.xml',weights_name{w}, segments{seg}, r), 'w');
                else
                    g = fopen(sprintf('xmls/Flexneri_%s_%s_rep%d.xml',weights_name{w}, segments{seg}, r), 'w');
                end                
                % keep track of the segment count for the nexus file name
                segmentcount = 1; segmentcount2=1;
                while ~feof(f_inf)
                    line = fgets(f_inf);
                    if ~isempty(strfind(line, 'insert_taxa'))
                    elseif contains(line, 'insert_sequence')
                        if startsWith(segs_files(seg).name,'LN')
                            fasta = fastaread(['data/Aln_Trees/LN624486_MDR_plasmid_aln/' segs_files(seg).name]);
                        elseif dataset_nr==1
                            fasta = fastaread(['data/Aln_Trees/Ss046_Sonnei_ref/' segs_files(seg).name]);
                        else
                            fasta = fastaread(['data/Aln_Trees/Flexneri_ref/' segs_files(seg).name]);
                        end

                        seq_length(seg) = length(fasta(1).Sequence);

                        for j = 1 : length(use_seqs)
                            ind = has_plasmid(use_seqs(j),seg);
                            if ind~=-1
                                fprintf(g, '\t\t<sequence id="%s.%s" taxon="%s" totalcount="4" value="%s"/>\n',...
                                     segments{seg}, fasta(ind).Header,...
                                     fasta(ind).Header, fasta(ind).Sequence);
                            end
                        end
                    elseif ~isempty(strfind(line, 'insert_sampling_times'))
                        indices = use_seqs(has_plasmid(use_seqs,seg)>-1);
                        stimes = 'rem';

                        for j = 1 : length(indices)
                            stimes = [stimes ',' ids{indices(j)} '=' times{indices(j)} '-01'];
                        end 
                        fprintf(g, strrep(line, 'insert_sampling_times', strrep(stimes, 'rem,','')));
                    elseif contains(line, 'insert_weights')
                        fprintf(g, strrep(line, 'insert_weights', num2str(slen)));
                    elseif ~isempty(strfind(line, 'initial_Ne'))
                        fprintf(g, strrep(line, 'initial_Ne', num2str(exprnd(1))));
                    elseif contains(line, 'insert_nexus_name')
                        line = strrep(line, 'insert_nexus_name', sprintf('sim_%d.%s.alignment.nexus', i,segment_names{seg}) );
                        fprintf(g, line);
                    else
                        fprintf(g, '%s', line);
                    end
                end

                fclose(f_inf); fclose(g);
            end
        end
    end
end

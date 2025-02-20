% makes simulation xmls for the coalescent with reassortment
clear

rng(156456454);

% define the number of repetitions
nr_reps = 50;

% rebuild the xml dirs
system('rm -r xmls');
system('mkdir xmls');

% define the number of samples
nr_samples = 100;

% define the sampling interval
sampling_interval = 5;

% define the rate shift points
rate_shifts = [200];


% file that keeps track of the ne and reassortment rates
h = fopen('rates.csv', 'w'); 
% fprintf(h, 'run,%splasmidTransfer\n', sprintf('popSize.%d,',1:length(rate_shifts)+1));
fprintf(h, 'run,popSize,plasmidTransfer\n');

% define params of the lognormal distribution of the reassortment rate
mean_rea = 0.01;
sigma_rea = 0.5;
mu_rea = log(mean_rea) - sigma_rea^2/2;

slen=[8000 200 100 50];                


segment_names = {'core1','plasmid1','plasmid2', 'plasmid3'};

% make nr reps number of xmls
for i = 1 : nr_reps
    f_sim = fopen('sim_template.xml');
    
    % sample the Ne and reassortment rates
    
    Ne = normrnd(3,1);
%     for j = 1 : length(rate_shifts)
%         Ne(j+1) = normrnd(Ne(j),1);
%     end
%     plot(Ne); hold on
    Ne = exprnd(100,1);

    reassortment = lognrnd(mu_rea, sigma_rea);
    fprintf(h, '%d,%s%.12f\n', i, sprintf('%.12f,',Ne), reassortment);
    
    % open the simulation xml
    g = fopen(sprintf('xmls/sim_%d.xml', i), 'w');
    
    use_sample = false(nr_samples, length(segment_names));
    
    while ~feof(f_sim)
        line = fgets(f_sim);
        if contains(line, 'insert_taxa')
            fprintf(g,'\t\t\t<taxonSet idref="taxonSet.%s"/>\n',segment_names{1});
            for j=2:length(segment_names)
                fprintf(g,'\t\t\t<plasmidTaxonSet idref="taxonSet.%s"/>\n',segment_names{j});
            end
        elseif contains(line, 'insert_trees') 
            time = zeros(nr_samples,1);
            for j = 1 : nr_samples
                time(j) = rand*sampling_interval;
            end

            
            for j=1:length(segment_names)
                fprintf(g,'\t\t\t<tree id="%s.tree" spec="Tree">\n',segment_names{j});
                fprintf(g,'\t\t\t\t<taxonset spec="TaxonSet" id="taxonSet.%s">\n',segment_names{j});
                if j==1
                    for k = 1 : floor(nr_samples)
                        use_sample(k,j) = true;
                        fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', k);
                    end
                else
                    for k = 1 : nr_samples
                        if rand < 1/j
                            use_sample(k,j) = true;
                            fprintf(g,'\t\t\t\t<taxon spec="Taxon" idref="t%d"/>\n', k);
                        end
                    end  
                end
                fprintf(g,'\t\t\t</taxonset>\n');

                fprintf(g,'\t\t\t\t<trait spec="TraitSet" traitname="date-backward" id="traitSet.%s">\n', segment_names{j});
                indices = find(use_sample(:,j));
                for k = 1:length(indices)
                    if k==length(indices)
                        fprintf(g,'\t\t\t\t\tt%d=%f\n', indices(k), time(indices(k)));
                    else
                        fprintf(g,'\t\t\t\t\tt%d=%f,\n', indices(k), time(indices(k)));
                    end
                end
                fprintf(g,'\t\t\t\t\t<taxa idref="taxonSet.%s"/>\n', segment_names{j});
                fprintf(g,'\t\t\t\t</trait>\n');
                fprintf(g,'\t\t\t</tree>\n');
            end



        elseif contains(line, 'insert_sampling_times')
            
            fprintf(g,'\t\t\t<traitSet idref="traitSet.%s"/>\n', segment_names{1});
        elseif contains(line, 'insert_seq_sim')
            for j=1:length(segment_names)
                fprintf(g,'\t\t\t<simulationObject id="align.%s" spec="SimulatedAlignment" outputFileName="$(filebase).%s.alignment.nexus" sequenceLength="%d" tree="@%s.tree">\n', segment_names{j}, segment_names{j},slen(j), segment_names{j});
                fprintf(g,'\t\t\t\t<siteModel spec="SiteModel">\n');
                fprintf(g,'\t\t\t\t\t<mutationRate spec="RealParameter" value="0.0005"/>\n');
                fprintf(g,'\t\t\t\t\t<substModel spec="JukesCantor"/>\n');
                fprintf(g,'\t\t\t\t</siteModel>\n');
                fprintf(g,'\t\t\t</simulationObject>\n');
            end

        elseif contains(line, 'insert_logger')
            for j=1:length(segment_names)
                fprintf(g,'\t\t\t<logger spec="Logger" fileName="$(filebase).%s.tree" logEvery="1">\n', segment_names{j});
                fprintf(g,'\t\t\t\t<log idref="%s.tree"/>\n', segment_names{j});
                fprintf(g,'\t\t\t</logger>\n');
            end      
            
        elseif contains(line, 'insert_Ne')                        
            fprintf(g, strrep(line, 'insert_Ne', num2str(Ne)));
        elseif contains(line, 'insert_rateShifts')
            fprintf(g, strrep(line, 'insert_rateShifts', num2str(rate_shifts)));

        elseif contains(line, 'insert_reassortment')
            fprintf(g, strrep(line, 'insert_reassortment', num2str(reassortment)));
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f_sim); fclose(g);
    %% make 3 replicates of the  network inference xmls
    for r = 1 : 3
        % build the inference xml
        f_inf = fopen('inf_template.xml');


        % open the inference xml
        g = fopen(sprintf('xmls/inf_%d_rep%d.xml', i, r), 'w');

        % keep track of the segment count for the nexus file name
        segmentcount = 1;segmentcount2=1;
%                 for j = 1 : nr_samples
%                     fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
%                 end             

        while ~feof(f_inf)
            line = fgets(f_inf);
            if ~isempty(strfind(line, 'insert_taxa'))
            elseif ~isempty(strfind(line, 'insert_sampling_times'))
                for j = 1 : nr_samples
                    if j==nr_samples
                        fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                    else
                        fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                    end
                end 
            elseif contains(line, 'insert_weights')
%                 fprintf(g, strrep(line, 'insert_weights', num2str(slen)));
                fprintf(g, strrep(line, 'insert_weights', num2str([1 1 1 1])));
            elseif contains(line, 'insert_rateShifts')
                fprintf(g, strrep(line, 'insert_rateShifts', num2str(rate_shifts)));

            elseif ~isempty(strfind(line, 'initial_Ne'))
                fprintf(g, strrep(line, 'initial_Ne', num2str(exprnd(1))));
            elseif ~isempty(strfind(line, 'initial_plasmidTransferRate'))
                fprintf(g, strrep(line, 'initial_plasmidTransferRate', num2str(exprnd(1))));
            elseif ~isempty(strfind(line, 'insert_sim_file_name'))
                line = strrep(line, 'insert_sim_file_name', sprintf('sim_%d', i) );
                fprintf(g, line);
            else
                fprintf(g, '%s', line);
            end
        end

        fclose(f_inf); fclose(g);
    end
    
    
    %% make 3 replicates of the tree inference xmls
    for r = 1 : 3

        for seg = 1 :length(segment_names)
            % build the inference xml
            f_inf = fopen('individualtrees_template.xml');

            % open the inference xml
            g = fopen(sprintf('xmls/inf_%d_%s_rep%d.xml', i, segment_names{seg}, r), 'w');

            % keep track of the segment count for the nexus file name
            segmentcount = 1; segmentcount2=1;

            while ~feof(f_inf)
                line = fgets(f_inf);
                if ~isempty(strfind(line, 'insert_taxa'))
                elseif ~isempty(strfind(line, 'insert_sampling_times'))
                    seg_samples = find(use_sample(:,seg));
                    for j = 1 : length(seg_samples)
                        
                        if j==length(seg_samples)
                            fprintf(g,'\t\t\t\tt%d=%f\n', seg_samples(j), time(seg_samples(j)));
                        else
                            fprintf(g,'\t\t\t\tt%d=%f,\n', seg_samples(j), time(seg_samples(j)));
                        end
                    end 
                elseif contains(line, 'insert_rateShifts')
                    fprintf(g, strrep(line, 'insert_rateShifts', num2str(rate_shifts)));

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
fclose(h);
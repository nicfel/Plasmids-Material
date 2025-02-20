% makes simulation xmls for the coalescent with reassortment
clear

rng(1);

% define the number of repetitions
nr_reps = 100;

% rebuild the xml dirs
system('rm -r xmls');
system('mkdir xmls');

% define the number of samples
nr_samples = 100;

% define the sampling interval
sampling_interval = 20;

% file that keeps track of the ne and reassortment rates
h = fopen('rates.csv', 'w');fprintf(h, 'run,Ne,plasmidTransfer\n');

% define params of the lognormal distribution of the Ne
mean_ne = 5;
sigma_ne = 0.5;
mu_ne = log(mean_ne) - sigma_ne^2/2;

% define params of the lognormal distribution of the reassortment rate
mean_rea = 0.2;
sigma_rea = 0.5;
mu_rea = log(mean_rea) - sigma_rea^2/2;

% make nr reps number of xmls
for i = 1 : nr_reps
    f_sim = fopen('sim_template.xml');
    
    % sample the Ne and reassortment rates
    Ne = lognrnd(mu_ne,sigma_ne);
    reassortment = lognrnd(mu_rea, sigma_rea);
    fprintf(h, '%d,%.12f,%.12f\n', i, Ne, reassortment);
    
    % open the simulation xml
    g = fopen(sprintf('xmls/sim_%d.xml', i), 'w');
    
    while ~feof(f_sim)
        line = fgets(f_sim);
        if ~isempty(strfind(line, 'insert_taxa'))
            for j = 1 : nr_samples
                fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
            end             
        elseif ~isempty(strfind(line, 'insert_sampling_times'))
            time = zeros(nr_samples,1);
            for j = 1 : nr_samples
                time(j) = rand*sampling_interval;
                if j==nr_samples
                    fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                else
                    fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                end
            end 
        elseif ~isempty(strfind(line, 'insert_Ne'))
            fprintf(g, strrep(line, 'insert_Ne', num2str(Ne)));
        elseif ~isempty(strfind(line, 'insert_reassortment'))
            fprintf(g, strrep(line, 'insert_reassortment', num2str(reassortment)));
        else
            fprintf(g, '%s', line);
        end
    end
    fclose(f_sim); fclose(g);
    
    % make 3 replicates
    for r = 1 : 3
        % build the inference xml
        f_inf = fopen('inf_template.xml');


        % open the inference xml
        g = fopen(sprintf('xmls/inf_%d_rep%d.xml', i, r), 'w');

        % keep track of the segment count for the nexus file name
        segmentcount = 1;segmentcount2=1;

        while ~feof(f_inf)
            line = fgets(f_inf);
            if ~isempty(strfind(line, 'insert_taxa'))
                for j = 1 : nr_samples
                    fprintf(g,'\t\t\t\t<taxon spec="Taxon" id="t%d"/>\n', j);
                end             
            elseif ~isempty(strfind(line, 'insert_sampling_times'))
                for j = 1 : nr_samples
                    if j==nr_samples
                        fprintf(g,'\t\t\t\tt%d=%f\n', j, time(j));
                    else
                        fprintf(g,'\t\t\t\tt%d=%f,\n', j, time(j));
                    end
                end 
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

end
fclose(h);
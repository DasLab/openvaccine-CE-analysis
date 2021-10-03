 function output = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, norm_pk, output_dir, output_pdf, fit_type );
% output = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, norm_pk, output_dir, output_pdf );
%
% INPUTS
%   data_dir    = input directory -- must contain
%                      sample_nucleotide_filename.csv
%                      platenumber_filename.csv
%                      processed_data/ subdirectory
%                 as created by Do Soon Kim's scripts.
%   pk_nt_bounds = [2xN] array of bounds of possible peaks. Last one is
%                    assumed to be peak of interest. 
%   norm_pk     = which of N peaks to use as reference. Give 0 to estimate
%                   total mRNA (total area minus first N-1 reference peaks)
%                   for normalization.
%   output_dir  = output directory for PDF's, .csv's. (default 'output')
%   output_pdf  = output PDF's (takes time).
%   fit_type    = 'expfit' or 'expfit_varyamp'
%
% OUTPUT
%   output = struct with fields like halflife, kdeg, rel_error, etc.
%
% (C) R. Das, Stanford University, 2021

% how much to go below or above each peak to capture area; scale by w
% outputted from findpeaks()
PKWIDTH_BELOW = 1; PKWIDTH_ABOVE = 2;

MIN_NT_FOR_TOTAL_AREA = 10;
MAX_NT_FOR_TOTAL_AREA = 1000;


% nucleotide bounds for 25-mer reference, P4-P6 reference, nLuc mRNA -- 
if ~exist( 'pk_nt_bounds', 'var' ) | isempty( pk_nt_bounds); pk_nt_bounds = [20 30; 150 300; 820 1000]; end
if ~exist( 'norm_pk', 'var' ); norm_pk = 0; end
if ~exist( 'output_dir', 'var' ); output_dir = 'output'; end
if ~exist( output_dir, 'dir' ) mkdir( ['output/',output_dir] ); end;
if ~exist( 'output_pdf', 'var' ); output_pdf = 0; end
if ~exist( 'fit_type', 'var' ); fit_type = 'expfit'; end
NUM_PKS_CURATED = size( pk_nt_bounds,1);

[~,data_dir_name,~] = fileparts(data_dir);

% sample_nucleotide_filename.csv has info on all samples.
x = readtable([data_dir,'/sample_nucleotide_filename.csv'],'Delimiter',',')
dx = table2cell( x );
all_sample_names = dx(:,1);
all_nucleotide = dx(:,2);
timepoint = cell2mat(dx(:,3));
plate = cell2mat(dx(:,4));
filenumber = cell2mat(dx(:,5));

if length(unique(all_nucleotide))>1 % could have PSU, m5C, etc.
    % to make sample names, add nucleotide tag.
    for m = 1:length(all_sample_names)
        all_sample_names{m} = [all_sample_names{m},'-',all_nucleotide{m}];
    end
end

sample_names = unique(all_sample_names,'stable')
kdeg = NaN*ones(1,length(sample_names));
halflife = NaN*ones(1,length(sample_names));
rel_error = NaN*ones(1,length(sample_names));

% platenumber_filename matches "platenumber" (actually chip number) with
%    filename outputted from Bioanalyzer and Do Soon's scripts.
x = readtable([data_dir,'/platenumber_filename.csv'],'Delimiter',',')
dx = table2cell( x );
plate_number = cell2mat(dx(:,1));
plate_file_names = dx(:,2);

tau_fit=zeros(1,length(sample_names));

% create .csv of output
outfile = sprintf('output/%s_sample_nt_filename_QUANT_DATA.csv',output_dir);
fid = fopen( outfile, 'w' );
fprintf( fid, 'Sample,Nucleotide,Timepoint,Plate,FileNumber,FracIntactNorm,FracIntact,TotalArea,Total_mRNA');
for i = 1:NUM_PKS_CURATED; fprintf( fid, ',PeakArea%d',i); end;
fprintf( fid,'\n');

for m = 1:length( sample_names )
    
    % sample by sample
    sample_name = sample_names{m};
    sample_idx = find(strcmp(all_sample_names, sample_name) & ~isnan(plate) );

    % which .csv files have processed data for this sample?
    num_samples = length(sample_idx);
    degtimes= timepoint(sample_idx);
    csvfile = {};
    for n = 1:num_samples
        csvfiles{n} = '';        
        if isnan( plate(sample_idx(n)) ) continue;end;
        file_name = plate_file_names{ find( plate_number == plate(sample_idx(n))) };
        csvfiles{n} = sprintf('%s/processed_data/nts-%s_Sample%d.csv',data_dir,file_name,filenumber(sample_idx(n)));
        if ~exist( csvfiles{n}, 'file' )
            csvfiles{n} = sprintf('%s/Analysis_DSK/%s_%d_nts.csv',data_dir,file_name,filenumber(sample_idx(n)));
        end
    end
    %if isempty( csvfiles{n} ); continue; end;
    if (num_samples == 0 ) continue; end;
    
    %  Let's do it!
    set( figure(1), 'position', [86    57   473   898] );
    clf
    pk_present = zeros(1,NUM_PKS_CURATED);
    total_area = zeros(1,NUM_PKS_CURATED);
    loc_curated = zeros(NUM_PKS_CURATED,num_samples);
    w_curated = zeros(NUM_PKS_CURATED,num_samples);
    peak_area = zeros(NUM_PKS_CURATED,num_samples);
    backgd_area = zeros(NUM_PKS_CURATED,num_samples);
    peak_area_backsub = zeros(NUM_PKS_CURATED,num_samples);
    
    for n = 1:num_samples
        subplot(length(csvfiles),1,n);
        
        % format of .csv files
        x = readtable(csvfiles{n});
        dx = cell2mat(table2cell( x ));
        time = dx(:,2);
        fu   = dx(:,3);
        nts  = dx(:,4);
        
        % plot the data.
        plot( time, fu); hold on
        nt_marks = [10:10:30,40:20:100,150:50:200,300:100:1000, 1200:200:1600];
        time_marks = interp1( nts, time, nt_marks );
        fprintf('%s --> %d %d\n', sample_name, m, n );
        gp = find( ~isnan(time_marks) & (time_marks > [0 time_marks(1:end-1)]) );
        set(gca,'xtick',time_marks(gp),'xticklabel',nt_marks(gp),'xticklabelrot',90);
        xlim( [20 40] )
        if ( n == 1) h=title( sample_name ); set(h,'interp','none');end;

         % find prominent peaks
        [pks,loc,w,p] = findpeaks( fu, 'MinPeakProminence',0.5 );
        % look for peaks that match the N specified in pk_nt_bounds.
        for i = 1:NUM_PKS_CURATED
            peak_idx = find( nts(loc) > pk_nt_bounds(i,1) & nts(loc) < pk_nt_bounds(i,2) & w<30);
            if length( peak_idx ) > 0
                [~,max_idx] = max(p(peak_idx)); peak_idx = peak_idx(max_idx);
                %peak_idx = max(peak_idx);
                loc_curated(i,n) = loc(peak_idx);
                w_curated(i,n) = w(peak_idx);
            end
        end
        
        % some additional curation -- if we do not find the peak in first 
        %  sample, assume its not there. 
        % If the peak is present in first sample but we didn't find in
        %  current sample, use peak location/width in previous sample.
        for i = 1:NUM_PKS_CURATED
            if ( n == 1)
                if loc_curated(i,n)>0; pk_present(i) = 1; end;
            else
                if pk_present(i);
                    if loc_curated(i,n) == 0;
                        loc_curated(i,n) = loc_curated(i,n-1);
                        w_curated(i,n)   = w_curated(i,n-1);
                    end
                else
                    loc_curated(i,n) = 0;
                end;
            end
        end
        
        %  peak areas
        for q = 1:NUM_PKS_CURATED
            if ( loc_curated(q,n) == 0 ) continue; end;

            % following is somewhat arbitrary
            min_idx = floor(loc_curated(q,n) - PKWIDTH_BELOW *w_curated(q,n));
            max_idx = ceil( loc_curated(q,n) + PKWIDTH_ABOVE *w_curated(q,n));
            
            % plotting
            area( time([min_idx:max_idx]), fu([min_idx:max_idx]));
            area( [time(min_idx) time(max_idx)], [fu(min_idx) fu(max_idx)] );
            
            % straight up integration
            peak_area(q,n)   = sum(fu([min_idx:max_idx]));
            % trapezoid for background subtraction
            backgd_area(q,n) = (max_idx-min_idx) * (fu(min_idx)+fu(max_idx))/2;
            peak_area_backsub(q,n) = peak_area(q,n) - backgd_area(q,n);
        end

        % total area 
        min_idx = floor(interp1( nts, [1:length(fu)], MIN_NT_FOR_TOTAL_AREA ));
        if isnan( min_idx ) min_idx = 1; end;
        max_idx = ceil( interp1( nts, [1:length(fu)], MAX_NT_FOR_TOTAL_AREA ));
        total_area(n) = sum( fu([min_idx:max_idx]) );
    end
    
    total_mRNA = total_area-sum(peak_area_backsub(1:(NUM_PKS_CURATED-1),:));
    set(gcf, 'PaperPositionMode','auto','color','white');
    if output_pdf; export_fig( sprintf('output/%s/Traces_Sample%02d_%s.pdf',output_dir,m,sample_name)); end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    clf
    subplot(2,1,1)
    set(gcf, 'PaperPositionMode','auto','color','white');
    plot( degtimes, total_mRNA,'o-');
    hold on
    plot( degtimes,peak_area_backsub','o-');
    legend_titles{1} = 'total-mRNA';
    for n = 1:NUM_PKS_CURATED; legend_titles{1+n} = ['Peak ',num2str(n)]; end;
    legend( legend_titles );
    h = title( sprintf('%s\n%s',sample_name, data_dir_name) ); set(h,'interp','none')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Normalize data to prepare for fitting
    subplot(2,1,2);    
    if norm_pk > 0
        frac_intact = peak_area_backsub(3,:)./peak_area_backsub(norm_pk,:);
    elseif norm_pk == 0
        frac_intact = peak_area_backsub(3,:)./total_mRNA;
    elseif norm_pk < 0
        frac_intact = peak_area_backsub(3,:);
    end
    d_norm = frac_intact/frac_intact(1);
            
    NBOOTSTRAP = 100;
    which_timepts = [1:length(d_norm)]; % find( d_norm >-0.05);

    % Do the fit
    p0 = [1 1]; %fit_type = 'expfit_varyamp';
    p = fminsearch( fit_type, p0, [], degtimes(which_timepts), d_norm(which_timepts)' );
    tau_fit(m) = p(1); 
    amp_fit(m) = 1;
    if strcmp(fit_type,'expfit_varyamp'); amp_fit(m) = p(2); end;
    
    % bootstrap
    tau_boot = [];
    if NBOOTSTRAP > 0
        %fprintf( 'Calculating %d bootstraps...\n',NBOOTSTRAP);
        for n = 1:NBOOTSTRAP
            boot_set = randi( length(which_timepts), length(which_timepts), 1 );
            p = fminsearch( fit_type, p0, [], degtimes(which_timepts(boot_set)), d_norm(which_timepts(boot_set))' );
            tau_boot(n) = p(1);
        end
        tau_fit_mean_boot(m) = mean(tau_boot);
        tau_fit_err_boot(m) = std(tau_boot);
    end
    
    plot( degtimes(which_timepts), d_norm(which_timepts),'ko' ); hold on
    finetimes = [0:0.1:25];
    plot( finetimes, amp_fit(m) * exp(-finetimes/tau_fit(m)), 'k-' );
    plot( [0 25],[0 0],'color',[0.5 0.5 0.5]);
    xlabel( 'Degradation time (hours)');
    ylabel( 'Frac. intact' );
    xlim([0 25]);
    rel_error(m) = tau_fit_err_boot(m)/tau_fit(m);
    kdeg(m) = 1/tau_fit(m);
    halflife(m) = log(2)*tau_fit(m);
    title( sprintf('k_{deg} = %5.3f+/-%5.3f /h. Half-life = %5.3f+/-%5.3f h. Norm peak: %d',kdeg(m),kdeg(m)*rel_error(m),...
        halflife(m),halflife(m)*rel_error(m), norm_pk));
    set(gcf, 'PaperPositionMode','auto','color','white');
    if output_pdf; export_fig( sprintf('output/%s/ExpFit_Sample%02d_%s.pdf',output_dir,m,sample_name)); end;
    %pause;
    
    % Output to .csv file
    for n = 1:num_samples
        i = sample_idx(n);
        fprintf( fid, '%s,%s,%f,%d,%d,%8.5f,%8.5f,%8.5f,%8.5f',...
            all_sample_names{i},all_nucleotide{i},timepoint(i),plate(i),filenumber(i),d_norm(n),frac_intact(n),total_area(n),total_mRNA(n));
        for j = 1:NUM_PKS_CURATED; fprintf( fid, ',%8.5f',peak_area_backsub(j,n)); end;
        fprintf( fid,'\n');
    end
        
end

fprintf( 'Created... %s\n', outfile );
fclose(fid);

% Summary graph of half lives (bar plot)
set(figure(3),'position',[ 553    35   560   420]);
clf
barh( halflife ); set(gca,'ytick',1:length(sample_names),'yticklabel',sample_names,'ticklabelinterp','none');set(gca,'ydir','reverse')
hold on
for m = 1:length(halflife)
    plot( halflife(m)*[1-rel_error(m), 1+rel_error(m)], m*[1 1],'k','linew',2);
end
maxval = max(halflife(find(~isnan(halflife))));
xlim([0 maxval*1.05])
h = title( sprintf('Half life (h), %s',data_dir_name) ); set(h,'interp','none');
set(gcf, 'PaperPositionMode','auto','color','white');
if output_pdf; export_fig( sprintf('output/%s_AllHalfLife_BarPlot.pdf',output_dir)); end;

% output summary to .csv
outfile = sprintf('output/%s_exp_fits_MATLAB.csv',output_dir);
fid = fopen( outfile, 'w' );
fprintf(fid,'Samples,Nucleotide,k_deg per hour,k_deg_err per hour\n');
for m = 1:length( sample_names )    
    % sample by sample
    sample_name = sample_names{m};
    sample_idx = find(strcmp(all_sample_names, sample_name) );
    i = sample_idx(1);    
    fprintf(fid,'%s,%s,%8.5f,%8.5f\n',...
        sample_name,all_nucleotide{i},kdeg(m),kdeg(m)*rel_error(m));
end
fprintf( 'Created... %s\n', outfile );
fclose(fid);


% package useful variables into output struct
output = struct();
output.kdeg = kdeg;
output.halflife = halflife;
output.rel_error = rel_error;
output.data_dir_name = data_dir_name;
output.output_dir = output_dir;
output.sample_names = sample_names;

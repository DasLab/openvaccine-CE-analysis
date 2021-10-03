pk_nt_bounds = [20 30; 150 300; 820 1000]; 
output_pdf = 0;
%% 
%% 

%% 12-02 (4-timepoint survey), P4-P6, mRNA
data_dir = 'data/120220_Bioanalyzer';
pk_norm = 2; output_dir = 'output_12-02_P4P6norm'; 
output_1202_P4P6norm = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf );

%% 12-02 (4-timepoint survey), P4-P6, mRNA -- expftivaryamp
data_dir = 'data/120220_Bioanalyzer';
pk_norm = 2; output_dir = 'output_12-02_P4P6norm_expfitvaryamp';
output_1202_P4P6norm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf, 'expfit_varyamp');

%%  12-02 (4-timepoint survey), use total-mRNA
data_dir = 'data/120220_Bioanalyzer';
pk_norm = 0; output_dir = 'output_12-02_totalmRNAnorm';
output_1202_totalmRNAnorm = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf );

%%  12-02 (4-timepoint survey), use total-mRNA -- expfitvaryamp
data_dir = 'data/120220_Bioanalyzer';
pk_norm = 0; output_dir = 'output_12-02_totalmRNAnorm';
output_1202_totalmRNAnorm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf,'expfit_varyamp' );



%% 12-10 (10-timepoint survey, no P4P6 spikein), use total-mRNA
data_dir = 'data/121020_Bioanalyzer';
pk_norm = 0; output_dir = 'output_12-10_totalmRNAnorm';
output_1210_totalmRNAnorm = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf );

%% 12-10 (10-timepoint survey, no P4P6 spikein), use total-mRNA Try to vary amplitudes.
data_dir = 'data/121020_Bioanalyzer';
pk_norm = 0; output_dir = 'output_12-10_totalmRNAnorm_expfitvaryamp';
output_1210_totalmRNAnorm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds, pk_norm, output_dir, output_pdf,'expfit_varyamp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6 nLuc texts (natural RNA) (10-15)
output_pdf = 1;
data_dir = 'data/101520_mRNA_BioAnalyzer';
pk_norm = 2; output_dir = 'output_10-15_6NLuc_P4P6norm'; 
pk_nt_bounds_6NLuc = [20 30; 150 250; 700 1000]; 
output_1015_6NLuc_P4P6norm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_6NLuc, pk_norm, output_dir, output_pdf,'expfit_varyamp'  );

%%
output_pdf = 1;
data_dir = 'data/101520_mRNA_BioAnalyzer';
pk_norm = 0; output_dir = 'output_10-15_6NLuc_totalmRNAnorm'; 
pk_nt_bounds_6NLuc = [20 30; 150 250; 700 1000]; 
output_1015_6NLuc_totalmRNAnorm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_6NLuc, pk_norm, output_dir, output_pdf,'expfit_varyamp'  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% m5C, PSU tests (10-22)
output_pdf = 1;
data_dir = 'data/102120_PSU_5mC_BioAnalyzer';
pk_norm = 2; output_dir = 'output_10-22_6Nluc_PSU_m5C_P4P6norm'; 
pk_nt_bounds_6NLuc = [20 30; 150 250; 700 1000]; 
output_1022_modRNA_6NLuc_P4P6_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_6NLuc, pk_norm, output_dir, output_pdf,'expfit_varyamp' );

%%
output_pdf = 1;
data_dir = 'data/102120_PSU_5mC_BioAnalyzer';
pk_norm = 0; output_dir = 'output_10-22_6Nluc_PSU_m5C_totalmRNAnorm'; 
pk_nt_bounds_6NLuc = [20 30; 150 300; 700 1000]; 
output_1022_modRNA_6NLuc_totalmRNAnorm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_6NLuc, pk_norm, output_dir, output_pdf,'expfit_varyamp' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UTR combo tests (11-30)
output_pdf = 0;
data_dir = 'data/113020_Bioanalyzer';
pk_norm = 2; output_dir = 'output_11-30_varyUTR_P4P6norm'; 
pk_nt_bounds_varyUTR = [20 30; 250 300; 700 1400]; 
output_1130_varyUTR_P4P6_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_varyUTR, pk_norm, output_dir, output_pdf,'expfit_varyamp' );

%%
output_pdf = 1;
data_dir = 'data/113020_Bioanalyzer';
pk_norm = 0; output_dir = 'output_11-30_varyUTR_P4P6norm'; 
pk_nt_bounds_varyUTR = [20 30; 250 300; 700 1400]; 
output_1130_varyUTR_totalmRNAnorm_expfitvaryamp = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_varyUTR, pk_norm, output_dir, output_pdf,'expfit_varyamp' );

%%
output_pdf = 0;
data_dir = 'data/113020_Bioanalyzer';
pk_norm = 0; output_dir = 'output_11-30_varyUTR_P4P6norm_noexptfitvaramp'; 
pk_nt_bounds_varyUTR = [20 30; 250 300; 700 1400]; 
output_1130_varyUTR_totalmRNAnorm = analyze_bioanalyzer_data( data_dir, pk_nt_bounds_varyUTR, pk_norm, output_dir, output_pdf );


%%
compare_half_life(output_1202_P4P6norm,output_1202_totalmRNAnorm);
%%
compare_half_life(output_1202_P4P6norm,output_1202_totalmRNAnorm_expfitvaryamp);
%%
compare_half_life(output_1210_totalmRNAnorm,output_1210_totalmRNAnorm_expfitvaryamp);
%%
compare_half_life(output_1202_P4P6norm,output_1210_totalmRNAnorm_expfitvaryamp);

%%  m5C/PSU data --> P4P6 norm vs. total area norm
compare_half_life(output_1022_modRNA_6NLuc_P4P6_expfitvaryamp,output_1022_modRNA_6NLuc_totalmRNAnorm_expfitvaryamp);

%% compile all data on one bar plot
sample_names = output_1202_P4P6norm.sample_names;
sample_colors = {[0.5 0.5 0.5],[0.8 0.3 0.3],[0.5 0.2 0.8],[0.5 0.2 0.8],...
    [0.9 0.5 0.2], [0.9 0.5 0.2], [0.9 0.5 0.2], [0.9 0.5 0.2],...
    [0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],[0 0.5 0],...
    [0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8],[0.2 0.3 0.8]};

% Rationale -- included P4-P6 normalization in 12/02/20 data set, so go
% ahead and use it. For the 12/10/20 repeat with more timepoints, neglected
% to include P4-P6 but, for 12/02, exptfit (varyamp) matches well to
% P4-P6-normalized.
outputs =  {output_1202_P4P6norm,output_1210_totalmRNAnorm_expfitvaryamp};
[halflife_mean,halflife_err] = half_life_compilation_plot( sample_names, sample_colors, outputs );

%%
%outfile = sprintf('averaged_k_deg_MATLAB_2020only_PROPAGATE_BOOTSTRAP_ERRORS.csv');
outfile = sprintf('averaged_k_deg_MATLAB_2020only.csv');
fid = fopen( outfile, 'w' );
fprintf(fid,'Samples,Nucleotide,k_deg per hour,k_deg_err per hour\n');
nucleotide = 'PSU';
kdeg = log(2) ./ halflife_mean
rel_error = halflife_err./halflife_mean;
for m = 1:length( sample_names )    
    % sample by sample
    sample_name = sample_names{m};
    fprintf(fid,'%s,%s,%8.5f,%8.5f\n',...
        sample_name,nucleotide,kdeg(m),kdeg(m)*rel_error(m));
end
fprintf( 'Created... %s\n', outfile );
fclose(fid);


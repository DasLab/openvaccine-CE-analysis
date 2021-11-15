#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy import optimize
import math
import matplotlib.pyplot as plt
import seaborn as sns

def analyze_bioanalyzer_data(data_dir, output_filename, sample_map_file = 'sample_nucleotide_filename.csv', plate_file_name = 'platenumber_filename.csv', 
    pk_nt_bounds=np.array([[20,30],[150,300],[820,1000]]), norm_pk=0, output_dir='output/', output_pdf=False, fit_type='expfit', remove_outliers = True, full_spike = False):
    '''
    INPUTS
        data_dir = (string) input directory, must contain:
                sample_nucleotide_filename.csv
                platenumber_filename.csv
                processed_data/ subdirectory with nts-converted csvs
        output_filename = (string) output filename for fit and calculation .csv and .pdfs
        pk_nt_bounds = [2xN] array of bounds of possible peaks. Last one assumed to be peak of interest.
        norm_pk = (int) which of N peaks to use as reference. 0 estimates based on total mRNA (total area minus first N-1 reference peaks) for normalization
        output_dir = output directory for PDF, csv outfiles (default = 'output/')
        output_pdf = (booleans) outputs PDF files, takes time
        fit_type = 'expfit' or 'expfit_varyamp', expfit fits y = e^(bx) and _varyamp fits y=Ae^(bx)

    OUTPUT
        output = dataframe with fields like halflife, kdeg, rel_err, etc.

    (C) D.S. Kim and R. Das, Stanford University, 2021 '''
    
    if os.path.exists(output_dir):
        pass
    else:
        os.mkdir(output_dir)

    print('Analyzing samples for: '+str(data_dir)+'...')
    
    #reading in samples and the corresponding filenames
    samples_df = pd.read_csv(data_dir+sample_map_file).dropna(axis=0)
    plates_df = pd.read_csv(data_dir+plate_file_name).dropna(axis=0)
    filenames_dict = dict(zip(plates_df['Plate_Number'],plates_df['File_Name']))

    all_samples = samples_df.to_dict('records')
    print(all_samples)

    for i, sample in enumerate(tqdm(all_samples)):
        print(sample)
        # assigning .csv file that contains the processed data for sample
        plate_num = sample['Plate']
        sample_num = str(sample['FileNumber']).split('.')[0]
        sample['filename'] = 'nts-'+filenames_dict[plate_num]+'_Sample'+sample_num+'.csv'
        
        # now read in data per sample and add to dictionary
        temp_df = pd.read_csv(data_dir+'/processed_data/'+sample['filename'])
        sample['time'] = np.array(temp_df['Time'], dtype=float)
        sample['fu'] = np.array(temp_df['Value'], dtype=float)
        sample['nt'] = np.array(temp_df['Nucleotides'], dtype=float)

        # finding peaks, prominence and width are set arbitrarily
        sample['peak_idx'], properties = find_peaks(sample['fu'], prominence=1, width=0.1, distance=30)
        # print(plate_num)
        # print(sample_num)
        # print(sample['peak_idx'])
        sample['widths'] = properties['widths']

        # filter for only peaks that fall within bounds in pk_nt_bounds
        # adding only for CoV-2 samples here:
        # print(sample['Sample'].split('_')[0])
        # if (sample['Sample'].split('_')[0]=='CoV-2-TTG-TTGfull-dSL1-3'):
        #     pk_nt_bounds = np.array([[20,30],[230,300],[1300,1600]])
        #     print(sample['Sample'].split('_')[0])
        #     print('triggered new pk_nt_bounds: '+ str(pk_nt_bounds))
        num_pks = len(pk_nt_bounds)
        loc_curated = []
        w_curated = []
        for q, idx in enumerate(sample['peak_idx']):
            for r in np.arange(num_pks):
        #         test_sample['nt'][idx]
                # for spike protein, increased max width to 50
                if (full_spike == True):
                    max_width = 50
                else:
                    max_width = 30
                if (sample['nt'][idx]>pk_nt_bounds[r][0] and sample['nt'][idx]<pk_nt_bounds[r][1] and sample['widths'][q]<max_width):
                    loc_curated.append(idx)
                    w_curated.append(sample['widths'][q])

        # store the peak indicies (loc_curated), associated widths (w_curated), and the nucleotides where peaks occur (peak_nts)
        loc_curated = np.array(loc_curated)
        w_curated = np.array(w_curated)

        sample['loc_curated'] = np.array(loc_curated)
        sample['w_curated'] = np.array(w_curated)
        sample['peak_nts'] = np.array(sample['nt'][loc_curated])

        # assume that the last peak is the peak of interest--if not found, then use w, loc from before
        # and if the peak doesn't exist in the first sample, assume that it's not there.

        if (full_spike == False): ### for most very precise RNAs, keep things as is
            if sample['Timepoint']==0.0:
                pass
            elif sample['nt'][loc_curated[-1]]>pk_nt_bounds[-1][0] and sample['nt'][loc_curated[-1]]<pk_nt_bounds[-1][1]:
                pass
            else:
                loc_curated = all_samples[i-1]['loc_curated']
                w_curated = all_samples[i-1]['w_curated']

        if (full_spike == True): ### if full spike protein mRNA (and thus likely diffuse)
            if sample['Timepoint']==0.0:
                pass
            else:
                timepoint0_dict = next((sub for sub in all_samples if ((sub['Sample'] == sample['Sample']) & (sub['Nucleotide']==sample['Nucleotide']) &(sub['Timepoint']==0.0))), None)
                print(timepoint0_dict)
                loc_curated = timepoint0_dict['loc_curated']
                w_curated = timepoint0_dict['w_curated']

        sample['loc_curated'] = np.array(loc_curated)
        sample['w_curated'] = np.array(w_curated)
        sample['peak_nts'] = np.array(sample['nt'][loc_curated])

        #calculate peak areas
        #somewhat arbitary here, capturing the search space
        pkwidth_below = 1
        pkwidth_above = 1

        min_idx = [math.floor(x) for x in (loc_curated-pkwidth_below*w_curated)]
        max_idx = [math.ceil(x) for x in (loc_curated+pkwidth_above*w_curated)]

        peak_areas = []
        backgd_areas = []
        peak_areas_backsub = []
        for minidx, maxidx in zip(min_idx, max_idx):
            # straight up integration for peak area
            area = sum(sample['fu'][minidx:maxidx])
            # background sub done by trapezoid area calc
            backgd_area = (maxidx-minidx)*((sample['fu'][minidx]+sample['fu'][maxidx])/2)
            # subtracting background
            peak_bkgd_subtract = area-backgd_area
            # add to lists
            peak_areas.append(area)
            backgd_areas.append(backgd_area)
            peak_areas_backsub.append(peak_bkgd_subtract)

        # assign lists to each sample
        sample['peak_area'] = np.array(peak_areas)
        sample['backgd_area'] = np.array(backgd_areas)
        sample['peak_area_backsub'] = np.array(peak_areas_backsub)

        # calculate total area for lane
        # first interpolate for nt closest to minimum and max window for total area
        MIN_NT_FOR_TOTAL_AREA = 20
        if (full_spike==False):
            MAX_NT_FOR_TOTAL_AREA = 1000
        if (full_spike==True):
            MAX_NT_FOR_TOTAL_AREA = 4600

        # deviates slightly from rhiju--he uses fu, I use nt in the interpolation. 
        # those two arrays should be equal length so no real change
        min_idx = math.floor(interp1d(sample['nt'], np.arange(len(sample['nt'])))(MIN_NT_FOR_TOTAL_AREA))
        max_idx = math.floor(interp1d(sample['nt'], np.arange(len(sample['nt'])))(MAX_NT_FOR_TOTAL_AREA))

        if np.isnan(min_idx):
            min_idx = 1

        total_area = sum(sample['fu'][min_idx:max_idx])
        sample['total_area'] = total_area

    # output to one df
    all_samples_df = pd.DataFrame(all_samples)
    all_samples_df.to_csv('checkpoint.csv')
    sample_df_list = []
    sample_fit_list = []

    for sample_nuc_pair in pd.unique(tuple(zip(all_samples_df['Sample'], all_samples_df['Nucleotide']))):
        sample_fit_dict = {}
        
        sample = sample_nuc_pair[0]
        sample_nucleotide = sample_nuc_pair[1]
        sample_fit_dict['Samples'] = sample
        sample_fit_dict['Nucleotide'] = sample_nucleotide

        print(sample)
        
        temp_df = all_samples_df[(all_samples_df['Sample']==sample) & (all_samples_df['Nucleotide']==sample_nucleotide)]
        peak_area_backsub = np.stack(temp_df['peak_area_backsub'].to_numpy())
        total_mRNA = temp_df['total_area'].to_numpy()

        # normalization depends on user input. 0: normalize by total, <0: no norm, >0 indicates normalizing peak
        if norm_pk == 0:
            frac_intact = np.array(peak_area_backsub[:,-1]/total_mRNA)
        elif norm_pk < 0:
            frac_intact = np.array(peak_area_backsub[:,-1])
        elif norm_pk > 0:
            frac_intact = np.array(peak_area_backsub[:,-1]/peak_area_backsub[:,norm_pk-1])
        frac_norm = (frac_intact/frac_intact[0]).flatten()

        temp_df.loc[:,'frac_intact'] = frac_intact
        temp_df.loc[:,'frac_norm'] = frac_norm

        #append results thus far to the comprehensive df
        sample_df_list.append(temp_df)

        # doing an exponential fit on the samples, both with and without option of varying amplitude
        expfit_f = lambda p,t,y: (abs(p[0]*np.exp(-1*t/p[1])-y)).sum()

        times = temp_df['Timepoint'].to_numpy()

        # do the fit!
        p0 = np.array([1,1])
        p = optimize.fmin(func=expfit_f, x0=p0, args = (times, frac_norm), disp=False)
        tau_fit = p[1]

        if fit_type == 'expfit_varyamp':
            amp = p[0]
        elif fit_type == 'expfit':
            amp = 1

        # bootstrap to calculate error
        N_BOOTSTRAP = 1000
        tau_boot = []
        for i in np.arange(N_BOOTSTRAP):
            bootstrap_inds = np.random.choice(a=len(times), size=len(times))
            fit_t = times[bootstrap_inds]
            fit_frac = frac_norm[bootstrap_inds]
            p = optimize.fmin(func=expfit_f, x0=p0, args=(fit_t,fit_frac), disp=False)
            tau_boot.append(p[1])
        tau_fit_mean_boot = np.mean(tau_boot)
        tau_fit_err_boot = np.std(tau_boot)

        if (remove_outliers == True):
            # let's filter for outliers, where cutoff is determined as a Z score of 3
            tau_boot = np.array(tau_boot)
            tau_fit_filtered_boot = tau_boot[np.where(abs(tau_boot-tau_fit_mean_boot)<=3*tau_fit_err_boot)]
            tau_fit_mean_boot = np.mean(tau_fit_filtered_boot)
            tau_fit_err_boot = np.std(tau_fit_filtered_boot)

        # from the fits, calculate kdeg and half-life
        rel_error = tau_fit_err_boot/tau_fit
        kdeg = 1/tau_fit
        halflife = np.log(2)*tau_fit

        # add to fit dictionary
        sample_fit_dict['k_deg per hour'] = kdeg
        sample_fit_dict['k_deg_err per hour'] = kdeg*rel_error
        sample_fit_dict['halflife (hours)'] = halflife
        sample_fit_dict['halflife_err (hours)'] = halflife*rel_error
        sample_fit_dict['amplitude'] = amp
        if (remove_outliers == True):
            sample_fit_dict['tau_boot'] = tau_fit_filtered_boot
        else:
            sample_fit_dict['tau_boot'] = tau_boot
        # add dictionary to list
        sample_fit_list.append(sample_fit_dict)

    compiled_df = pd.concat(sample_df_list).reset_index()
    fit_summary_df = pd.DataFrame(sample_fit_list).reset_index()
    
   # script for output to csv and json
    compiled_df.to_csv(output_dir+output_filename+'_PYTHON_sample_nt_filename_QUANT_DATA.csv')
    compiled_df.to_json(output_dir+output_filename+'_PYTHON_sample_nt_filename_QUANT_DATA.json')

    # also saving summary of fits to csv
    fit_summary_df.to_csv(output_dir+output_filename+'_exp_fits_PYTHON.csv')
    fit_summary_df.to_json(output_dir+output_filename+'_exp_fits_PYTHON.json')


    # outputting plots!
    if (output_pdf==True):
        for sample, nucleotide in zip(fit_summary_df['Samples'], fit_summary_df['Nucleotide']):    
            plt.clf()
            plt.figure(figsize=(8,12))

            ########################
            ### FIRST subPLOT: peak areas over time points
            ########################

            plt.subplot(2,1,1)

            raw_data = compiled_df[(compiled_df['Sample']==sample) & (compiled_df['Nucleotide']==nucleotide)]
            timepoints = np.array(raw_data['Timepoint'])
            total_areas = np.array(raw_data['total_area'])
            peak_areas = np.column_stack(raw_data['peak_area_backsub'])

            sns.scatterplot(x=timepoints, y=total_areas, marker='o', s=100, legend=False, linewidth=2, 
                edgecolor='blue', facecolor='none')
            plt.plot(timepoints, total_areas, label='Total mRNA', linewidth=1, color='blue')
            
            custom_colors = ['#7b113a', '#150e56', '#1597bb', '#8fd6e1']
            for i, peak_area in enumerate(peak_areas):
                sns.scatterplot(x=timepoints, y=peak_area, marker='o', s=100, legend=False, linewidth=1, edgecolor=custom_colors[i], 
                    facecolor='none')
                plt.plot(timepoints, peak_area, label='Peak '+str(i+1), linewidth=1, color=custom_colors[i])
                
            plt.xlabel('Time (hours)', fontsize=14)
            plt.ylabel('mRNA amount (area)', fontsize=14)
            plt.title('%s \n %s' %(sample, nucleotide),fontsize=12)
            plt.xticks(fontsize=13)
            plt.yticks(fontsize=13)
            plt.legend()

            ########################
            ### SECOND subPLOT: exponential fits for half life calculation
            ########################
            # grab kdeg and amplitude for exponential function fit

            plt.subplot(2,1,2)
            
            select_df = fit_summary_df[(fit_summary_df['Samples']==sample) & (fit_summary_df['Nucleotide']==nucleotide)]
            kdeg = float(select_df['k_deg per hour'])
            kdeg_err = float(select_df['k_deg_err per hour'])
            halflife = float(select_df['halflife (hours)'])
            halflife_err = float(select_df['halflife_err (hours)'])
            amplitude = float(select_df['amplitude'])
            
            # and now grab the fraction intact values and the corresponding time points
            frac_intacts = np.array(raw_data['frac_norm'])
            
            # let's plot it!
            # plt.figure(figsize=(8,12))
            sns.scatterplot(x=timepoints, y=frac_intacts, marker='o', s=100, legend=False, linewidth=1, edgecolor='black', facecolor='none')
            plt.xlabel('Time (hours)', fontsize=14)
            plt.ylabel('Fraction Intact', fontsize=14)
            plt.title('kdeg = %.3f +/- %.3f /h. Half-life = %.3f +/- %.3f h. Norm peak: %i' 
                  %(kdeg, kdeg_err, halflife, halflife_err, norm_pk),
                 fontsize=14)
            plt.xticks(fontsize=13)
            plt.yticks(fontsize=13)
            
            times = np.arange(0,max(timepoints), 0.05)
            fit_vals = amplitude*np.exp(-1*kdeg*times)
            plt.plot(times, fit_vals, '--' ,linewidth = 2, color='black')

            ### export and save figure
            plt.tight_layout()
            plt.savefig(output_dir+'ExpFit_'+'NormPeak'+str(norm_pk)+fit_type+'_'+sample+'_'+nucleotide+'.pdf')

            ########################
            ### THIRD PLOT: exponential fits for half life calculation
            ########################
            plt.clf()
            plt.figure(figsize=(10,20))
            plot_suptitle = plt.suptitle(sample+'_'+nucleotide, x=0.5, y=1.015, fontsize=14, fontweight='bold')
            
            for i, row in enumerate(raw_data.itertuples()):
                
                plt.subplot(len(raw_data), 1, i+1)
                
                nts = np.array(row.nt)
                fus = np.array(row.fu)
                times = np.array(row.time)
                
                
                # just for the spike protein analysis (4000 nucleotides)
                if (full_spike == True):
                    nt_marks = np.array([10,20,30,40,60,80,100,150,200,300,400,500,1000,2000,3000,4000,4500])

                    times_curated = times[100:int(len(times))]
                    nts_curated = nts[100:int(len(nts))]
                    fus_curated = fus[100:int(len(fus))]


                else:
                    nt_marks = np.array([10,20,30,40,60,80,100,150,200,300,400,500,600,700,800,900,1000,1200,1400])
                    times_curated = times[100:int(len(times)/2.2)]
                    nts_curated = nts[100:int(len(nts)/2.2)]
                    fus_curated = fus[100:int(len(fus)/2.2)]

                time_marks = interp1d(nts, times, fill_value="extrapolate")(nt_marks)

                # figure out the fill_between indices, in times
                # first figure out min and max idx, same as above
                pkwidth_below = 1; pkwidth_above = 2;
                
                min_idx = [math.floor(x) for x in (np.array(row.loc_curated)-pkwidth_below*np.array(row.w_curated))]
                max_idx = [math.ceil(x) for x in (np.array(row.loc_curated)+pkwidth_above*np.array(row.w_curated))]
                print(sample)
                print(min_idx)
                print(max_idx)

                plt.plot(times_curated, fus_curated)
                
                #color in the peaks, rotate color
                custom_colors = ["#ffc93c","#9ddfd3","#31326f"]
                num_peaks = np.arange(len(pk_nt_bounds))
                for min_i, max_i, ccolor, count in zip(min_idx, max_idx, custom_colors, num_peaks):
                    plt.fill_between(times[min_i:max_i], fus[min_i:max_i], color=ccolor, label='Peak '+str(count+1))
                    if count+1==len(num_peaks):
                        plt.fill_between(x=times[min_i:max_i], y1 = np.linspace(fus[min_i], fus[max_i], len(times[min_i:max_i])), color='gray', label='Background')
                    else:
                        plt.fill_between(x=times[min_i:max_i], y1 = np.linspace(fus[min_i], fus[max_i], len(times[min_i:max_i])), color='gray')
                    plt.legend()
                plt.xticks(ticks=time_marks, labels=nt_marks, fontsize=13)
                plt.yticks(fontsize=13)
                plt.ylabel("Signal (a.u.)", fontsize=14)
                if (full_spike==True):
                    plt.xlim(22,54)
                plt.ylim(0,25)
                plt.title("Timepoint: %s hrs" %(row.Timepoint), fontsize=14)

            plt.xlabel("Nucleotides", fontsize=15)
            plt.tight_layout()
            plt.savefig(output_dir+'Traces_'+sample+'_'+nucleotide+'.pdf', bbox_inches='tight', bbox_extra_artists=[plot_suptitle])

    print('Completed!')
    return compiled_df, fit_summary_df
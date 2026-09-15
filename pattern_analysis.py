import os
import pandas as pd
import numpy as np
from re import search
from ssa_schlicke_25 import do_full_fit_of_pattern_notime


# directories where the data can be found and the figures should be saved
path_data = 'C:/Users/dummy/insilico_biopsies_from_that_weird_paper'  # 'C:/Users/vroni/Documents/Uni/Promotion/Papers/Paper_1/github_code/insilico_biopsies_from_that_weird_paper'
path_figures = 'figures'

# input ts for model-PSD and 2PCF that should be used
fixed_ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# true parameters (order: gamma, D, K, maxT)
V = [0.3, 0.005, 1, 6]  # minimum proliferation, minimum diffusion
W = [0.7, 0.15, 1, 6]  # max prolif, max diff
X = [0.7, 0.005, 1, 6]  # max prolif, min diff
Y = [0.3, 0.15, 1, 6]  # min prolif, max diff
parameters = {'V': V, 'W': W, 'X': X, 'Y': Y}


# for nice output of warnings
class bcolors:
    WARNING = '\033[93m'
    ENDC = '\033[0m'


if __name__ == "__main__":
    print('Starting SSA Optimization.')

    for t_fix in fixed_ts:
        print(f"{bcolors.WARNING}Starting Optimization for t = {t_fix}...{bcolors.ENDC}")

        # file where the estimates should be saved
        save_popts_file = f'{path_data}/RESULTS/popts_t={t_fix}.npz'

        # pre-allocate dictionaries to save the results
        gammas = {}
        Ds = {}
        ts = {}
        r_squared_2PCF_dict = {}
        r_squared_PSD_dict = {}
        gamma_err_dict = {}
        D_err_dict = {}
        t_err_dict = {}
        N_dict = {}
        area_dict = {}
        results = pd.DataFrame()
        if os.path.isdir(path_data):
            print(f'Parameter estimates will be saved in {save_popts_file}')
            if not os.path.exists(os.path.dirname(save_popts_file)):
                os.makedirs(os.path.dirname(save_popts_file))
            for file in os.listdir(path_data):
                # ignore files that start with '.' - potential MAC issue
                if file.startswith('.'):
                    continue
                if file.endswith('.csv'):
                    print(file)
                    # obtain group name from file name

                    match = search(r'([A-Z]{1,2})', file)
                    if match:
                        group = match.group(1)
                    else:
                        continue

                    # obtain the true parameters for this group
                    true_gamma = parameters[group][0]
                    true_D = parameters[group][1]
                    true_t = parameters[group][3]
                    # perform the fit using the current fixed t
                    gamma_fit, D_fit, r_squared_2PCF, r_squared_PSD, gamma_err, D_err, N, area = do_full_fit_of_pattern_notime(
                        path_data, file, path_figures, [true_gamma, true_D, true_t], 0, 0,
                        fix_t=t_fix, whole_domain=True, cut_N=False, cut_2PCF=False)

                    ts[group] = t_fix
                    # save the fit results
                    gammas[group] = gamma_fit
                    Ds[group] = D_fit
                    r_squared_2PCF_dict[group] = r_squared_2PCF
                    r_squared_PSD_dict[group] = r_squared_PSD
                    gamma_err_dict[group] = gamma_err
                    D_err_dict[group] = D_err
                    N_dict[group] = N
                    area_dict[group] = area

            # check if there already exists a dictionary whose values will be overwritten, and values from new groups will be added
            if os.path.isfile(save_popts_file):
                print(f'The file {save_popts_file} already exists and will be extended.')
                npzfile = np.load(save_popts_file, allow_pickle=True)
                gammas = npzfile['gammas'].item() | gammas
                Ds = npzfile['Ds'].item() | Ds
                ts = npzfile['ts'].item() | ts
                r_squared_2PCF_dict = npzfile['r_squared_2PCF'].item() | r_squared_2PCF_dict
                r_squared_PSD_dict = npzfile['r_squared_PSD'].item() | r_squared_PSD_dict
                gamma_err_dict = npzfile['gamma_err'].item() | gamma_err_dict
                D_err_dict = npzfile['D_err'].item() | D_err_dict
                t_err_dict = npzfile['t_err'].item() | t_err_dict
                N_dict = npzfile['N'].item() | N_dict
                area_dict = npzfile['area'].item() | area_dict

            np.savez(save_popts_file, gammas=gammas, Ds=Ds, ts=ts, r_squared_2PCF=r_squared_2PCF_dict,
                     r_squared_PSD=r_squared_PSD_dict, gamma_err=gamma_err_dict, D_err=D_err_dict,
                     t_err=t_err_dict, N=N_dict, area=area_dict)

            print('Optimization finished.')

# experiment type can be selected by commenting the lines after "experiment i", i = 1, 2 in the code

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm


# settings
path_data = 'C:/Users/dummy/insilico_biopsies_from_that_weird_paper'

key_exceptions = []  # optional, if some parameter combinations should be left out

gamma_fit = 'pwlinear'  # 'quadra'  'pwlinear'

thr_rsquared = 0.9  # SSA quality of fit standard

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


def sort_data(pathname, popt_filenames, thr_rsquared, key_exceptions):
    # dictionaries for the results  -- pre-allocation for automatic import
    gamma_estimates = {}
    D_estimates = {}
    r_squared_PSD = {}
    r_squared_2PCF = {}
    gamma_errs = {}
    D_errs = {}
    keylist = []
    for file in popt_filenames:
        npzfile = np.load(pathname + file, allow_pickle=True)  #os.path.join(path_matlab, file))
        opt_gammas = npzfile['gammas'].item()
        opt_Ds = npzfile['Ds'].item()
        saved_gamma_errs = npzfile['gamma_err'].item()
        saved_D_errs = npzfile['D_err'].item()
        saved_r_squared_PSD = npzfile['r_squared_PSD'].item()
        saved_r_squared_2PCF = npzfile['r_squared_2PCF'].item()
        for key in opt_gammas:
            if key in key_exceptions:
                pass
            elif key in gamma_estimates:
                gamma_estimates[key].append(opt_gammas[key])
            else:
                keylist.append(key)
                gamma_estimates[key] = [opt_gammas[key]]
        for key in opt_Ds:
            if key in key_exceptions:
                pass
            elif key in D_estimates:
                D_estimates[key].append(opt_Ds[key])
            else:
                D_estimates[key] = [opt_Ds[key]]
        for key in saved_gamma_errs:
            if key in key_exceptions:
                pass
            elif key in gamma_errs:
                gamma_errs[key].append(saved_gamma_errs[key])
            else:
                gamma_errs[key] = [saved_gamma_errs[key]]
        for key in saved_D_errs:
            if key in key_exceptions:
                pass
            elif key in D_errs:
                D_errs[key].append(saved_D_errs[key])
            else:
                D_errs[key] = [saved_D_errs[key]]
        for key in saved_r_squared_PSD:
            if key in key_exceptions:
                pass
            elif key in r_squared_PSD:
                r_squared_PSD[key].append(saved_r_squared_PSD[key])
            else:
                r_squared_PSD[key] = [saved_r_squared_PSD[key]]
        for key in saved_r_squared_2PCF:
            if key in key_exceptions:
                pass
            elif key in r_squared_2PCF:
                r_squared_2PCF[key].append(saved_r_squared_2PCF[key])
            else:
                r_squared_2PCF[key] = [saved_r_squared_2PCF[key]]

    # Sort x-values to ensure monotonic order
    parameters_subgroup = {k: parameters[k] for k in keylist if
                           k in parameters}  # select the subsection from GroupDictionary.parameters that contains the desired letter combinations
    gammas_sorted = sorted(parameters_subgroup, key=lambda k: parameters_subgroup[k][0])
    Ds_sorted = sorted(parameters_subgroup, key=lambda k: parameters_subgroup[k][1])
    D_over_gamma_true = {key: values[1] / values[0] for key, values in parameters_subgroup.items()}
    D_over_gamma_true_sorted = sorted(D_over_gamma_true, key=D_over_gamma_true.get)
    x1_values = []
    x2_values = []
    x3_values = []
    Ds_sorted_by_gammas = []
    Ds_sorted_by_Dovergamma = []
    for trial_set in range(len(labels)):
        x1_values.append([parameters_subgroup[k][0] for k in gammas_sorted])  # set of x-values for the gammas
        x2_values.append([parameters_subgroup[k][1] for k in Ds_sorted])  # set of x-values for the Ds
        x3_values.append(
            [parameters_subgroup[k][1] / parameters_subgroup[k][0] for k in
             D_over_gamma_true_sorted])  # set of x-values for D/gamma
        Ds_sorted_by_gammas.append(
            [parameters_subgroup[k][1] for k in gammas_sorted])  # true Ds in the same order as the true gammas
        Ds_sorted_by_Dovergamma.append(
            [parameters_subgroup[k][1] for k in
             D_over_gamma_true_sorted])  # true Ds in the same order as the true D/gamma

    # Extract y-, D/gamma-, R^2-values and errors in the same order as x1, x2 (ACHTUNG: y1_values und y2_values sind nicht in der selben Reihenfolge!!)
    y1_values = [gamma_estimates[k] for k in gammas_sorted]
    y2_values = [D_estimates[k] for k in Ds_sorted]
    r_squared_PSD_sortedbygamma = [r_squared_PSD[k] for k in gammas_sorted]
    r_squared_PSD_sortedbyD = [r_squared_PSD[k] for k in Ds_sorted]
    r_squared_2PCF_sortedbygamma = [r_squared_2PCF[k] for k in gammas_sorted]
    r_squared_2PCF_sortedbyD = [r_squared_2PCF[k] for k in Ds_sorted]
    gamma_errs_sortedbygamma = [gamma_errs[k] for k in gammas_sorted]
    D_errs_sortedbyD = [D_errs[k] for k in Ds_sorted]
    D_over_gamma_fit = []
    fitted_Ds_sorted_by_gammas = []
    for trial_set in range(len(labels)):
        D_over_gamma_fit.append(
            [D_estimates[k][trial_set] / (gamma_estimates[k][trial_set] + 1e-16) for k in D_over_gamma_true_sorted])
        fitted_Ds_sorted_by_gammas.append([D_estimates[k][trial_set] for k in
                                           gammas_sorted])  # fitted Ds in the same order as y1_values (fitted gammas by true gamma)

    # re-order the elements because the original order is too confusing
    y1_values = list(map(list, zip(*y1_values)))
    y2_values = list(map(list, zip(*y2_values)))
    r_squared_PSD_sortedbygamma = list(map(list, zip(*r_squared_PSD_sortedbygamma)))
    r_squared_PSD_sortedbyD = list(map(list, zip(*r_squared_PSD_sortedbyD)))
    r_squared_2PCF_sortedbygamma = list(map(list, zip(*r_squared_2PCF_sortedbygamma)))
    r_squared_2PCF_sortedbyD = list(map(list, zip(*r_squared_2PCF_sortedbyD)))
    gamma_errs_sortedbygamma = list(map(list, zip(*gamma_errs_sortedbygamma)))
    D_errs_sortedbyD = list(map(list, zip(*D_errs_sortedbyD)))
    # fitted_Ds_sorted_by_gammas = list(map(list, zip(*fitted_Ds_sorted_by_gammas)))

    if len(popt_filenames) == 1:
        pass

    # exclude zeros (zero means that the estimation was not possible) & fits with R^2 below a threshold
    disqualified_counter = np.zeros(len(labels), dtype=int)
    parameters_sorted_by_gamma = sorted(parameters_subgroup.items(), key=lambda item: item[1][0])

    # prepare lists to save R²s and fitting-stds
    mean_rsq_array = []
    rel_D_err_array = []
    rel_gamma_err_array = []

    for trial_set in range(len(labels)):
        del_idx1 = []  # gamma-related deletions
        del_idx2 = []  # D-related deletions
        for i in range(len(y1_values[trial_set])):
            # for each parameter: if the R² of PSD OR 2PCF is > thr_rsquared, it will be accepted (if a fit is available)
            if y1_values[trial_set][i] == 0 or (r_squared_PSD_sortedbygamma[trial_set][i] < thr_rsquared and
                                                r_squared_2PCF_sortedbygamma[trial_set][i] < thr_rsquared):
                del_idx1.append(i)
                ith_key = parameters_sorted_by_gamma[i][0]
                print('In the trial labelled "' + labels[trial_set] + '", a parameter fit had both R² < ' +
                      str(thr_rsquared) + ':  PSD-R² = ' +
                      str(round(r_squared_PSD_sortedbygamma[trial_set][i], 3)) + ', 2PCF-R² = ' +
                      str(round(r_squared_2PCF_sortedbygamma[trial_set][i], 3)) +
                      '. This occured for parameter combination ' +
                      f"{bcolors.WARNING}{ith_key}{bcolors.ENDC}" + '. Involved parameters: gamma = ' +
                      str(parameters_subgroup[ith_key][0]) + ', D = ' + str(parameters_subgroup[ith_key][1]))
                disqualified_counter[trial_set] += 1
            if y2_values[trial_set][i] == 0 or (
                    r_squared_PSD_sortedbyD[trial_set][i] < thr_rsquared and r_squared_2PCF_sortedbyD[trial_set][
                i] < thr_rsquared):  # r_squared_PSD_sortedbyD[trial_set][i] < thr_rsquared:
                del_idx2.append(i)
                # no need for a console output because if R² is too low for gamma it automatically also is for D

        y1_values[trial_set] = [y1_values[trial_set][j] for j in range(len(y1_values[trial_set])) if j not in del_idx1]
        x1_values[trial_set] = [x1_values[trial_set][j] for j in range(len(x1_values[trial_set])) if j not in del_idx1]
        gamma_errs_sortedbygamma[trial_set] = [gamma_errs_sortedbygamma[trial_set][j] for j in
                                               range(len(gamma_errs_sortedbygamma[trial_set])) if j not in del_idx1]
        r_squared_PSD_sortedbygamma[trial_set] = [r_squared_PSD_sortedbygamma[trial_set][j] for j in
                                                  range(len(r_squared_PSD_sortedbygamma[trial_set])) if
                                                  j not in del_idx1]
        r_squared_2PCF_sortedbygamma[trial_set] = [r_squared_2PCF_sortedbygamma[trial_set][j] for j in
                                                   range(len(r_squared_2PCF_sortedbygamma[trial_set])) if
                                                   j not in del_idx1]
        Ds_sorted_by_gammas[trial_set] = [Ds_sorted_by_gammas[trial_set][j] for j in
                                          range(len(Ds_sorted_by_gammas[trial_set])) if j not in del_idx1]

        y2_values[trial_set] = [y2_values[trial_set][j] for j in range(len(y2_values[trial_set])) if j not in del_idx2]
        x2_values[trial_set] = [x2_values[trial_set][j] for j in range(len(x2_values[trial_set])) if j not in del_idx2]
        D_errs_sortedbyD[trial_set] = [D_errs_sortedbyD[trial_set][j] for j in
                                       range(len(D_errs_sortedbyD[trial_set])) if j not in del_idx2]
        r_squared_PSD_sortedbyD[trial_set] = [r_squared_PSD_sortedbyD[trial_set][j] for j in
                                              range(len(r_squared_PSD_sortedbyD[trial_set])) if j not in del_idx2]
        r_squared_2PCF_sortedbyD[trial_set] = [r_squared_2PCF_sortedbyD[trial_set][j] for j in
                                               range(len(r_squared_2PCF_sortedbyD[trial_set])) if j not in del_idx2]

        # usually, if an entry in y1_values is 0, it is also 0 in y2_values, so this SHOULD work -- no guarantees
        D_over_gamma_fit[trial_set] = [D_over_gamma_fit[trial_set][j] for j in
                                       range(len(D_over_gamma_fit[trial_set])) if j not in del_idx1]
        x3_values[trial_set] = [x3_values[trial_set][j] for j in range(len(x3_values[trial_set])) if j not in del_idx1]
        Ds_sorted_by_Dovergamma[trial_set] = [Ds_sorted_by_Dovergamma[trial_set][j] for j in
                                              range(len(Ds_sorted_by_Dovergamma[trial_set])) if j not in del_idx1]

        # compute and print mean R²s of the remaining parameter combinations
        mean_rsquared = np.mean(r_squared_PSD_sortedbyD[trial_set])  # 0.5 * (np.mean(r_squared_PSD_sortedbyD[trial_set]) + np.mean(r_squared_2PCF_sortedbyD[trial_set]))
        mean_gamma_err = np.mean([gamma_errs_sortedbygamma[trial_set][i] / y1_values[trial_set][i] for i in range(len(y1_values[trial_set]))])
        mean_D_err = np.mean([D_errs_sortedbyD[trial_set][i] / y2_values[trial_set][i] for i in range(len(y2_values[trial_set]))])
        # print(f'In the accepted parameter combinations of trial "{labels[trial_set]}", the mean of PSD- and 2PCF-R² is '
        #       f'{mean_rsquared:.3f}. The mean fitting-std. during PSD- and 2PCF-fitting is given by D: '
        #       f'{mean_D_err:.3f}, gamma: {mean_gamma_err:.3f}')
        mean_rsq_array.append(r_squared_PSD_sortedbyD[trial_set])
        rel_gamma_err_array.append([gamma_errs_sortedbygamma[trial_set][i] / y1_values[trial_set][i] for i in range(len(y1_values[trial_set]))])
        rel_D_err_array.append([D_errs_sortedbyD[trial_set][i] / y2_values[trial_set][i] for i in range(len(y2_values[trial_set]))])

    # print mean R²s etc. for all trials
    print(f'\nOverall R² (PSD): {np.mean(flatten(mean_rsq_array)):.3f}, overall D-std: '
          f'{np.mean(flatten(rel_D_err_array)):.3f}, overall gamma-std: {np.mean(flatten(rel_gamma_err_array)):.3f}')

    return x1_values, y1_values, fitted_Ds_sorted_by_gammas, Ds_sorted_by_gammas, x2_values, y2_values  # xn: true parameter values, yn: estimated parameter values


def take_sample(fixed_ts, sample_size, x1_values, y1_values, fitted_Ds_sorted_by_gammas, true_Ds_sorted_by_gammas,
                x2_values, y2_values):
    # Der Output ist keine Liste von listen, sondern eine einzelne liste -> macht jetzt nicht mehr so viel Sinn streng nach t zu trennen
    N = sum(len(sublist) for sublist in x1_values)  # total number of value pairs (D, gamma)
    # np.random.seed(0)
    indices = np.random.randint(0, N, size=sample_size)
    x1_values_f = flatten(x1_values)
    y1_values_f = flatten(y1_values)
    fitted_Ds_sorted_by_gammas_f = flatten(fitted_Ds_sorted_by_gammas)
    true_Ds_sorted_by_gammas_f = flatten(true_Ds_sorted_by_gammas)
    x2_values_f = flatten(x2_values)
    y2_values_f = flatten(y2_values)
    x1 = []
    y1 = []
    fDsg = []
    tDsg = []
    x2 = []
    y2 = []
    ts = []
    # get the true and estimated gamma & D values at the specified locations
    for idx in indices:
        x1.append(x1_values_f[idx])
        y1.append(y1_values_f[idx])
        D_hat = fitted_Ds_sorted_by_gammas_f[idx]
        fDsg.append(D_hat)
        D = true_Ds_sorted_by_gammas_f[idx]
        tDsg.append(D)
        # _1 und _2 sind ja unterschiedlich sortiert, also muss ich den wahren D-Wert anders als über idx herausfinden, z.B. über eine Wert -> Index Suche
        idx_D = y2_values_f.index(D_hat)
        x2.append(x2_values_f[idx_D])
        y2.append(y2_values_f[idx_D])
        # get the points in time
        idx_copy = idx
        t_nr = 0
        for sublist in x1_values:
            if idx_copy < len(sublist):
                value = sublist[idx_copy]  # man braucht den Wert nicht mehr, ist vllt gut zum Debuggen
                break
            else:
                idx_copy -= len(sublist)
                t_nr += 1
        ts.append(fixed_ts[t_nr])
    return [x1], [y1], [fDsg], [tDsg], [x2], [y2], ts


# inverse quadratic form for the D-values
def inverse_quadratic(y, a, b, c):
    y = np.array(y)
    return a * np.sqrt(y) + b * y + c


def grouping_by_D(ascending_conv_Ds, fitted_gammas_sorted_by_D, x1_sorted_by_D):
    # aufteilen in Gruppen wo das convertierte D nur p * 100 % voneinander entfernt ist (rel. zur Gesamtrange)
    p = 0.1  # v18: 0.10  # 0.25  # 0.2
    # Ausreißer finden und entfernen
    ascending_conv_Ds, fitted_gammas_sorted_by_D, x1_sorted_by_D = find_and_eliminate_outliers(ascending_conv_Ds, fitted_gammas_sorted_by_D, x1_sorted_by_D)
    conv_D_range = max(ascending_conv_Ds) - min(ascending_conv_Ds)  # TODO: robuste Möglichkeit, Ausreißer auszuschließen, finden
    slice_idx = [0]  # um die indices an denen getrennt werden soll zu speichern
    D0 = ascending_conv_Ds[0]  # initialize D0
    for i in range(len(ascending_conv_Ds)):
        if ascending_conv_Ds[i] - D0 > p * conv_D_range:
            D0 = ascending_conv_Ds[i]
            slice_idx.append(i)

    slice_idx.append(len(ascending_conv_Ds))

    nr_groups = len(slice_idx) - 1  # bei der Anzahl der Gruppen muss man den ersten Index ignorieren

    conv_D_groups = []
    fitted_gammas_groups = []
    x1_sorted_groups = []

    # print('Groups (converted D-values): ')
    for i in range(len(slice_idx) - 1):
        conv_D_groups.append(ascending_conv_Ds[slice_idx[i]:slice_idx[i + 1]])
        fitted_gammas_groups.append(fitted_gammas_sorted_by_D[slice_idx[i]:slice_idx[i + 1]])
        x1_sorted_groups.append(x1_sorted_by_D[slice_idx[i]:slice_idx[i + 1]])

        # print(ascending_conv_Ds[slice_idx[i]:slice_idx[i + 1]])
        # if i == len(slice_idx) - 2:
        #     print('\n')

    return slice_idx, nr_groups, conv_D_groups, fitted_gammas_groups, x1_sorted_groups


def grouping_by_gamma(gamma_estimate, gamma_true, D_converted):
    # Jans Version von "grouping_by_D"
    gamma_estimate_sorted, gamma_true_sorted, D_converted_sorted = zip(
        *sorted(zip(gamma_estimate, gamma_true, D_converted)))

    slice_idx = [0]
    for i in range(1, len(gamma_estimate_sorted)):
        if (gamma_estimate_sorted[i] - gamma_estimate_sorted[i - 1]) > .06:
            slice_idx.append(i)
    slice_idx.append(len(gamma_estimate_sorted))

    gamma_true_groups = []
    gamma_estimate_groups = []
    D_converted_groups = []

    nr_groups = 0
    for i in range(len(slice_idx) - 1):
        tmp1 = []
        tmp2 = []
        tmp3 = []
        for j in range(slice_idx[i], slice_idx[i + 1]):
            tmp1.append(gamma_true_sorted[j])
            tmp2.append(gamma_estimate_sorted[j])
            tmp3.append(D_converted_sorted[j])
        gamma_true_groups.append(tmp1)
        gamma_estimate_groups.append(tmp2)
        D_converted_groups.append(tmp3)
        nr_groups += 1

    return slice_idx, nr_groups, D_converted_groups, gamma_estimate_groups, gamma_true_groups


def find_and_eliminate_outliers(ascending_conv_Ds, fitted_gammas_sorted_by_D, x1_sorted_by_D):
    """Findet Ausreißer in ascending_conv_Ds und passt alle Listen entsprechend an."""
    q1 = np.quantile(ascending_conv_Ds, 0.25)  # 1st Quartile
    q3 = np.quantile(ascending_conv_Ds, 0.75)  # 3rd Quartile
    iqr = np.sqrt((q3 - q1) ** 2)  # interquartile range
    outlier_lb = q1 - 1.5 * iqr  # lower bound of outliers
    outlier_ub = q3 + 1.5 * iqr  # upper bound of outliers
    ascending_conv_Ds_no_outliers = []
    fitted_gammas_sorted_by_D_no_outliers = []
    x1_sorted_by_D_no_outliers = []
    for idx, conv_D in enumerate(ascending_conv_Ds):
        if outlier_lb < conv_D < outlier_ub:
            ascending_conv_Ds_no_outliers.append(conv_D)
            fitted_gammas_sorted_by_D_no_outliers.append(fitted_gammas_sorted_by_D[idx])
            x1_sorted_by_D_no_outliers.append(x1_sorted_by_D[idx])
        else:
            print(f"{bcolors.WARNING}Outlier found and eliminated in converted Ds: {conv_D:.4f}{bcolors.ENDC}")
    return ascending_conv_Ds_no_outliers, fitted_gammas_sorted_by_D_no_outliers, x1_sorted_by_D_no_outliers


# find the hyperplane
def findHyperplane(X, Y, Z):
    # returns a,b,c in z=ax+by+c
    # assumes that len(X)=len(Y)=len(Z)
    A = np.array([[X[i], Y[i], 1] for i in range(len(X))])
    b = np.array(Z)
    res = np.linalg.solve(A.T @ A, A.T @ b)
    return res


# use the hyperplane
def linear_2vars(x, a, b, c):
    try:
        gamma = np.array(x[:, 0])
        D = np.array(x[:, 1])
    except:
        gamma = np.array(x[0])
        D = np.array(x[1])
    return a * gamma + b * D + c


def flatten(xss):
    return [x for xs in xss for x in xs]


def quadratic(x, m, n, o):
    x = np.array(x)
    return m * x ** 2 + n * x + o


def quadratic_fit(average_D_groups, popts_groups):
    popts_mno = np.zeros((3, 3))
    for param in range(3):
        popt_param, pcov, infodict, _, _ = curve_fit(quadratic, average_D_groups, popts_groups[:, param],
                                                     full_output=True)
        popts_mno[param, :] = popt_param
    return popts_mno


def piecewise_linear(x, all_x, all_y, const_ends=True):
    """
    :param x: die unabhängige Variable an der die Fkt. ausgewertet werden soll (array).
    :param all_x: array der x-Werte aller Bruchstellen.
    :param all_y: array der y-Werte aller Bruchstellen, sortiert wie all_x.
    :param const_ends: wenn 'True', werden die Bereiche [-infty, all_x[0]] und [all_x[-1], infty] gleich dem ersten/letzten all_y-Wert gesetzt.
    :return: Funktionswerte bei x.
    """
    if const_ends:
        incs = [0]  # für x < all_x[0] soll die Fkt. gleich all_y[0] sein
        ints = [all_y[0]]
        # add \infty to the list of breaking points, weil [-\infty, all_x[0]] und [all_x[-1], \infty] jetzt eigene Fkt.-Vorschriften haben
        breaks = [-np.inf] + all_x + [np.inf]
    else:
        incs = []
        ints = []
        # die Steigung der ersten und letzten Teilgeraden soll ins Unendliche beibehalten werden
        breaks = [-np.inf] + all_x[1:-1] + [np.inf]
    for i in range(len(all_x) - 1):
        # piecewise linear functions of the form y = incline * x + intersect
        incline = (all_y[i + 1] - all_y[i]) / (all_x[i + 1] - all_x[i])
        incs.append(incline)
        intersect = all_y[i] - incline * all_x[i]
        ints.append(intersect)
    if const_ends:
        incs.append(0)  # für x > all_x[-1] soll die Fkt. gleich all_y[-1] sein
        ints.append(all_y[-1])
    # conditions based on breakpoints
    conds = [(x >= breaks[i]) & (x < breaks[i + 1]) for i in range(len(breaks) - 1)]
    # Function list: each lambda captures the corresponding m and t
    funcs = [lambda z, m=incs[i], t=ints[i]: m * z + t for i in range(len(incs))]
    return np.piecewise(x, conds, funcs)


def D_over_gamma_curve(x, a, b, c, d, e, f, g):
    # x = [[D, gamma], ..., [D, gamma]]
    try:
        D = np.array(x[:, 0])
        gamma = np.array(x[:, 1])
    except:
        D = np.array(x[0])
        gamma = np.array(x[1])
    return (a * np.sqrt(D) + b * D + c) / (d * gamma + e * np.sqrt(D) + f * D + g)


def compute_RRMSE(rmse, mean_true_vals):
    """
    :param rmse: root mean squared error sqrt( 1/n * sum_{i=1}^n (y_i - hat{y}_i)² )
    :param mean_true_vals: Durchschnitt der wahren Werte 1/n * sum_{i=1}^n y_i
    :return: relative RMSE
    """
    return rmse / mean_true_vals


def compute_R_square(ssr, sq_ges):
    """
    :param ssr: sum of squared residuals sum_i (hat{y}_i - y_i)²
    :param sq_ges: Quadratsumme der wahren Werte sum_i (y_i - mean(y))²
    :return: R²
    """
    return 1 - ssr / sq_ges


if __name__ == "__main__":
    # experiment 1: if multiple fixed ts should be tested
    popt_files = []
    labels = []
    fixed_ts = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # [6]
    for t_fix in fixed_ts:
        popt_files.append(f'/RESULTS/popts_t={t_fix}.npz')
        labels.append(f't = {t_fix}')

    # arrays for the results
    popts_D_array = []
    rrmse_D_array = []
    rsq_D_array = []
    popts_gamma_array = []
    popts_gamma_quadr_array = []  # die m, n, o Parameter zu a, b, c falls quadratisch gefitted wird (nicht piecewise linear)
    rrmse_gamma_array = []
    rsq_gamma_array = []
    rrmse_gamma_reconversion_array = []
    rsq_gamma_reconversion_array = []

    # go through fitting results and get data in sorted manner
    x1_values, y1_values, fitted_Ds_sorted_by_gammas, true_Ds_sorted_by_gammas, x2_values, y2_values = sort_data(
        path_data, popt_files, thr_rsquared, key_exceptions)

    # experiment 2: take a sample from the data ignoring t-cohorts
    # sample_size = 30
    # x1_values, y1_values, fitted_Ds_sorted_by_gammas, true_Ds_sorted_by_gammas, x2_values, y2_values, ts = take_sample(fixed_ts, sample_size, x1_values, y1_values, fitted_Ds_sorted_by_gammas, true_Ds_sorted_by_gammas, x2_values, y2_values)
    # # wenn man sampelt muss man labels neu definieren
    # labels = [f'selection of {sample_size} estimates']

    # plots with fit_param on the x-axis, and true_param on the y ######################################################
    print(f'\nConversion function parameters for f: v_fit -> v_true')
    fig, axes = plt.subplots(1, 3)  # , figsize=(24, 8))
    fig.suptitle(f'In-Silico Data: Fitting results vs. true parameter values, R² > {thr_rsquared}')

    # retrieve & plot D
    # D: D_fit on x-axis, D_true on y
    markers = ['o', 's', '^', 'v', '<', '>', 'd', 'p', 'h', '*', 'x', '+']
    for trial_set in range(len(labels)):  # Number of lines
        general_xD_inv = np.linspace(0, max(y2_values[trial_set]), 100)
        axes[0].scatter(y2_values[trial_set], x2_values[trial_set], marker=markers[trial_set], alpha=0.5,
                        label=labels[trial_set])  # color='#808080',
        popt_D, pcov, infodict, _, _ = curve_fit(inverse_quadratic, y2_values[trial_set], x2_values[trial_set],
                                                 full_output=True)
        squared_errors = [i ** 2 for i in infodict['fvec']]
        rmse = np.sqrt(np.mean(squared_errors))
        rrmse = compute_RRMSE(rmse, np.mean(x2_values[trial_set]))
        rsq = compute_R_square(np.sum(squared_errors),
                               np.sum([(x2_values[trial_set][i] - np.mean(x2_values[trial_set])) ** 2 for i in
                                       range(len(x2_values[trial_set]))]))
        print(f'TRIAL {labels[trial_set]}: Fit results for D: a = {popt_D[0]}, b = {popt_D[1]}, c = {popt_D[2]}.')
        print(f'Root mean squared error: {rmse}')
        print(f'Relative root mean squared error: {rrmse}')
        print(f'R² = {rsq}\n')
        popts_D_array.append(
            [popt_D[0], popt_D[1], popt_D[2]])  # muss man so komisch machen weil curve_fit ein "array-Objekt" ausgibt
        rrmse_D_array.append(rrmse)
        rsq_D_array.append(rsq)
        axes[0].plot(general_xD_inv, inverse_quadratic(general_xD_inv, popt_D[0], popt_D[1], popt_D[2]))
    axes[0].set_ylabel("true D")
    axes[0].set_xlabel("estimated D")
    axes[0].set_title("D (diffusion)")
    axes[0].legend()
    axes[0].grid(True)

    # retrieve gamma

    # list to save results
    D_gamma_estimated_all_trials = []
    D_gamma_converted_all_trials = []
    D_gamma_true_all_trials = []

    for trial_set in range(len(labels)):
        # headline for the outlier reports and error outputs
        print(f'TRIAL {labels[trial_set]}: gamma: Conversion errors')

        # use popt_D from before to convert the Ds that are sorted by gamma
        converted_Ds = inverse_quadratic(fitted_Ds_sorted_by_gammas[trial_set], *popts_D_array[trial_set])

        # converted_Ds gruppieren in g Gruppen
        ascending_conv_Ds, fitted_gammas_sorted_by_D, x1_sorted_by_D, true_Ds_sorted_by_asc_c_D, fitted_Ds_sorted_by_asc_c_D = map(
            list, zip(*sorted(
                zip(list(converted_Ds), y1_values[trial_set], x1_values[trial_set], true_Ds_sorted_by_gammas[trial_set],
                    fitted_Ds_sorted_by_gammas[trial_set]))))

        # Nach konvertiertem D sortieren ###############################################################################
        slice_idx, nr_groups, conv_D_groups, fitted_gammas_groups, x1_sorted_groups = grouping_by_D(ascending_conv_Ds,
                                                                                                    fitted_gammas_sorted_by_D,
                                                                                                    x1_sorted_by_D)

        # Alternative: nach fitted gamma sortieren (Jan) ###############################################################
        # slice_idx, nr_groups, conv_D_groups, fitted_gammas_groups, x1_sorted_groups = grouping_by_gamma(fitted_gammas_sorted_by_D, x1_sorted_by_D, ascending_conv_Ds)

        # saving the popts and errors for the different groups
        popts_groups = []
        rmse_groups = []
        rrmse_groups = []
        rsq_groups = []
        group_sizes = []
        average_D_groups = []  # für's berechnen der Parameter-Beziehungen a(D), b(D), c(D) später
        squared_errors_groups = []

        for group in range(nr_groups):
            # curve_fit
            x_data = np.column_stack((fitted_gammas_groups[group], conv_D_groups[group]))
            try:
                popt_gamma, pcov, infodict, _, _ = curve_fit(linear_2vars, x_data, x1_sorted_groups[group], full_output=True, bounds=([0, -np.inf, -np.inf], [np.inf, np.inf, np.inf]))  # p0=[20, 50, 30])  #, # TODO: bei der erweiterten v18 brauchts das l.b.
                # popt_gamma, pcov, infodict, _, _ = curve_fit(alt_linear_2vars, x_data, x1_sorted_groups[group], full_output=True)
                squared_errors = [i ** 2 for i in infodict['fvec']]
            except RuntimeError:  # gedacht für den Fall "RuntimeError: Optimal parameters not found: The maximum number of function evaluations is exceeded."
                # einfach a, b, c von der letzten Gruppe wiederverwenden (popt_gamma, muss man nicht neu definieren)
                # squared errors neu ausrechnen
                fits = linear_2vars(x_data, popt_gamma[0], popt_gamma[1], popt_gamma[2])
                if len(fits) > 1:
                    squared_errors = [(fits[i] - x1_sorted_groups[group][i]) ** 2 for i in range(len(x1_sorted_groups[group]))]
                else:
                    squared_errors = (fits - x1_sorted_groups[group]) ** 2

            squared_errors_groups.append(squared_errors)
            rmse = np.sqrt(np.mean(squared_errors))
            rrmse = compute_RRMSE(rmse, np.mean(x1_sorted_groups[group]))  # TODO: sqrt nur im Zähler, kein ^2 im Nenner -> eq. 60 in https://www.sciencedirect.com/science/article/pii/S1364032115013258?via=ihub bzw https://www.marinedatascience.co/blog/2019/01/07/normalizing-the-rmse/
            rsq = compute_R_square(np.sum(squared_errors),
                                   np.sum([(x1_sorted_groups[group][i] - np.mean(x1_sorted_groups[group])) ** 2 for i in
                                            range(len(x1_sorted_groups[group]))]))

            popts_groups.append(popt_gamma)
            average_D_groups.append(np.mean([min(conv_D_groups[group]), max(conv_D_groups[group])]))

            rmse_groups.append(rmse)
            rrmse_groups.append(rrmse)
            rsq_groups.append(rsq)
            group_sizes.append(len(fitted_gammas_groups[group]))

        # gesamten Fehler noch ausrechnen (Durchschnitt der Fehler aller Gruppen)
        all_sq_errs = flatten(squared_errors_groups)
        mean_rmse = np.sqrt(np.mean(all_sq_errs))
        flat_x1_sorted_groups = flatten(x1_sorted_groups)
        mean_rrmse = compute_RRMSE(mean_rmse, np.mean(flat_x1_sorted_groups))
        mean_rsq = compute_R_square(np.sum(all_sq_errs),
                                    np.sum([(flat_x1_sorted_groups[i] - np.mean(flat_x1_sorted_groups)) ** 2 for i in
                                            range(len(flat_x1_sorted_groups))]))
        print(f'Overall root mean squared error: {mean_rmse}')
        print(f'Overall relative root mean squared error: {mean_rrmse}')
        print(f'R² = {mean_rsq}')

        rrmse_gamma_array.append(mean_rrmse)
        rsq_gamma_array.append(mean_rsq)

        popts_groups = np.array(
            popts_groups)  # each column is one parameter (a, b, c), the rows are the popt for the respective group

        # a(f(D)) = popts_mno[0, 0] * f²(D) + popts_mno[0, 1] * f(D) + popts_mno[0, 2] usw.
        # piecewise linear
        if gamma_fit == 'pwlinear':
            a_D = [piecewise_linear(d, average_D_groups, popts_groups[:, 0], const_ends=True) for d in
                   ascending_conv_Ds]
            b_D = [piecewise_linear(d, average_D_groups, popts_groups[:, 1], const_ends=True) for d in
                   ascending_conv_Ds]
            c_D = [piecewise_linear(d, average_D_groups, popts_groups[:, 2], const_ends=True) for d in
                   ascending_conv_Ds]
        # alternative: quadratic
        else:
            popts_mno = quadratic_fit(average_D_groups, popts_groups)
            a_D = [popts_mno[0, 0] * d ** 2 + popts_mno[0, 1] * d + popts_mno[0, 2] for d in ascending_conv_Ds]
            b_D = [popts_mno[1, 0] * d ** 2 + popts_mno[1, 1] * d + popts_mno[1, 2] for d in ascending_conv_Ds]
            c_D = [popts_mno[2, 0] * d ** 2 + popts_mno[2, 1] * d + popts_mno[2, 2] for d in ascending_conv_Ds]
            popts_gamma_quadr_array.append(popts_mno)  # (flatten(popts_mno))

        popts_gamma_array.append([a_D, b_D, c_D])

        sc = axes[1].scatter(fitted_gammas_sorted_by_D, x1_sorted_by_D, marker=markers[trial_set], alpha=0.5,
                             label=labels[trial_set])
        col = sc.get_facecolor()
        # cbar = plt.colorbar(sc)  # Attach colorbar to scatter plot
        # cbar.set_label("converted D-value")
        axes[1].set_ylabel("true gamma")
        axes[1].set_xlabel("estimated gamma")
        axes[1].set_title("gamma (proliferation)")
        axes[1].set_ylim(0, 1)
        axes[1].legend()
        axes[1].grid(True)

        x_data_all_Ds = np.column_stack((fitted_gammas_sorted_by_D, ascending_conv_Ds))
        f_of_gamma_fit = []
        for idx in range(len(ascending_conv_Ds)):
            general_x = np.linspace(min(fitted_gammas_sorted_by_D) - 0.05, max(fitted_gammas_sorted_by_D) + 0.05, 100)
            x_data = np.column_stack((general_x, [ascending_conv_Ds[idx]] * 100))
            axes[1].plot(general_x, linear_2vars(x_data, a_D[idx], b_D[idx], c_D[idx]), alpha=0.6, color=col)

            fogf = linear_2vars(x_data_all_Ds[idx, :], a_D[idx], b_D[idx], c_D[idx])
            f_of_gamma_fit.append(fogf)
        # RMSE und RRMSE berechnen (bissl komplizierter weil die Kurve ja nur so indirekt gefittet wurde)
        squared_errors = [(x1_sorted_by_D[i] - f_of_gamma_fit[i]) ** 2 for i in range(len(x1_sorted_by_D))]
        rmse = np.sqrt(np.mean(squared_errors))
        rrmse = compute_RRMSE(rmse, np.mean(x1_sorted_by_D))
        rsq = compute_R_square(np.sum(squared_errors), np.sum([(x1_sorted_by_D[i] - np.mean(x1_sorted_by_D)) ** 2 for i in range(len(x1_sorted_by_D))]))
        print(f'gamma: Re-conversion errors')
        print(f'Root mean squared error: {rmse}')
        print(f'Relative root mean squared error: {rrmse}')
        print(f'R² = {rsq}\n')

        rrmse_gamma_reconversion_array.append(rrmse)
        rsq_gamma_reconversion_array.append(rsq)

        # estimated Ds und gammas in der gleichen Reihenfolge (auch gleiche Reihenfolge wie die trues)
        D_gamma_estimated_all_trials.append(np.column_stack((fitted_Ds_sorted_by_asc_c_D, fitted_gammas_sorted_by_D)))
        # converted Ds und gammas in der gleichen Reihenfolge (auch gleiche Reihenfolge wie die trues)
        D_gamma_converted_all_trials.append(np.column_stack((ascending_conv_Ds, f_of_gamma_fit)))
        # true Ds und gammas in der gleichen Reihenfolge
        D_gamma_true_all_trials.append(np.column_stack((true_Ds_sorted_by_asc_c_D, x1_sorted_by_D)))

    # plot with est. D/gamma on the x-axis, and true D/gamma on the y
    # flatten the color-determining true D values to get global range
    all_D_values = [row[0] for trial in D_gamma_true_all_trials for row in trial]

    # normalize over the global value range
    norm = mcolors.Normalize(vmin=min(all_D_values), vmax=max(all_D_values))
    cmap = cm.cool  # Choose any colormap
    for trial_set in range(len(labels)):
        temp = [row[0] / row[1] for row in D_gamma_true_all_trials[trial_set]]
        temp2 = [row[0] / row[1] for row in D_gamma_estimated_all_trials[trial_set]]
        # wie temp & temp2 aber als np-arrays
        true_ratios = D_gamma_true_all_trials[trial_set][:, 0] / D_gamma_true_all_trials[trial_set][:, 1]
        estimated_ratios = D_gamma_estimated_all_trials[trial_set][:, 0] / D_gamma_estimated_all_trials[trial_set][:, 1]
        true_D_values = [row[0] for row in D_gamma_true_all_trials[trial_set]]  # Used for color

        colors = cmap(norm(true_D_values))  # RGBA colors

        # Step 4: Plot with color depending on D_value
        axes[2].scatter(temp2, temp, marker=markers[trial_set], alpha=0.5,
                        label=labels[trial_set], c=colors)
        # axes[2].scatter(temp2, temp, marker=markers[trial_set], alpha=0.5, label=labels[trial_set]) # color irrelevant

        if len(labels) == 1:  # sampling-Situation, in der D/gamma direkt retrieved werden
            popt, pcov, infodict, _, _ = curve_fit(D_over_gamma_curve, D_gamma_estimated_all_trials[trial_set], temp,
                                                   full_output=True)
            squared_errors = [i ** 2 for i in infodict['fvec']]
            rmse = np.sqrt(np.mean(squared_errors))
            rrmse = compute_RRMSE(rmse, np.mean(temp))
            rsq = compute_R_square(np.sum(squared_errors), np.sum([(temp[i] - np.mean(temp)) ** 2 for i in range(len(temp))]))
            print(f'TRIAL {labels[trial_set]}: Fit results for D/gamma')
            print(
                f'a = {popt[0]}, b = {popt[1]}, c = {popt[2]}, d = {popt[3]}, e = {popt[4]}, f = {popt[5]}, g = {popt[6]}.')
            print(f'Root mean squared error: {rmse}')
            print(f'Relative root mean squared error: {rrmse}')
            print(f'R² = {rsq}\n')

            # sort temp2 ([D_hat, gamma_hat]) in ascending order
            sorted_indices = np.argsort(estimated_ratios)
            ordered_temp2 = estimated_ratios[sorted_indices]  # estimated D/gamma
            ordered_D_gamma = D_gamma_estimated_all_trials[trial_set][sorted_indices]  # estimated [D, gamma]

            # temp ([D, gamma], also die wahren Werte) in der gleichen Reihenfolge, braucht man für Parity Plot
            ordered_temp = true_ratios[sorted_indices]

            D_over_gamma_converted = D_over_gamma_curve(ordered_D_gamma, *popt)

            axes[2].plot(ordered_temp2,
                         D_over_gamma_converted)  # Plot hat Knicke weil benachbarte D_hat/gamma_hat nicht zwingend durch benachbarte D_hats bzw. gamma_hats entstanden sind

    # Step 5: Add a single colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for colorbar to work without a mappable object
    cbar = plt.colorbar(sm, ax=axes[2])
    cbar.set_label('true D')  # Label for colorbar

    axes[2].set_ylabel("true D/gamma")
    axes[2].set_xlabel("estimated D/gamma")
    axes[2].set_title("D/gamma")
    axes[2].legend()
    axes[2].grid(True)

    plt.show()

    # Parity-Plots (das sind die mit der Winkelhalbierenden)
    fig, axes = plt.subplots(1, 3)
    if len(labels) == 1:  # unterschiedliche Plot-Titel je nach Experiment; früher war das einfach immer f'In-Silico Data: Parity Plots for converted fitted vs. true parameter values, R² > {thr_rsquared}'
        fig.suptitle(f'Parity Plots for converted fitted vs. true parameter values, t-cohorts unknown')
    else:
        fig.suptitle(f'Parity Plots for converted fitted vs. true parameter values, t-cohorts known')

    # D & gamma in nem loop
    rrmses = [rrmse_D_array, rrmse_gamma_reconversion_array]  # rrmse_gamma_array
    rsqs = [rsq_D_array, rsq_gamma_reconversion_array]
    xy_labels = ['D', 'gamma']
    for idx in [0, 1]:
        for trial_set in range(len(labels)):
            conv_temp = [row[idx] for row in D_gamma_converted_all_trials[trial_set]]
            true_temp = [row[idx] for row in D_gamma_true_all_trials[trial_set]]
            axes[idx].scatter(conv_temp, true_temp, marker=markers[trial_set], alpha=0.5, label=labels[trial_set])

        # RRMSE
        mean_rrmse = np.mean(rrmses[idx])
        std_rrmse = np.std(rrmses[idx])
        axes[idx].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f} +/- {std_rrmse:.3f}', transform=axes[idx].transAxes)

        # R²
        mean_rsq = np.mean(rsqs[idx])
        std_rsq = np.std(rsqs[idx])
        axes[idx].text(0.2, 0.05, f'R² = {mean_rsq:.3f} +/- {std_rsq:.3f}', transform=axes[idx].transAxes)

        if len(labels) == 1:  # mean & std machen keinen Sinn bei len(labels) == 1
            axes[idx].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f}', transform=axes[idx].transAxes)
            axes[idx].text(0.2, 0.05, f'R² = {mean_rsq:.3f}', transform=axes[idx].transAxes)
        else:
            axes[idx].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f} +/- {std_rrmse:.3f}', transform=axes[idx].transAxes)
            axes[idx].text(0.2, 0.05, f'R² = {mean_rsq:.3f} +/- {std_rsq:.3f}', transform=axes[idx].transAxes)

        # draw the bisection line (first, identify domain & range, then plot)
        axes[idx].set_xlim(left=0)
        axes[idx].set_ylim(bottom=0)
        xlim = axes[idx].get_xlim()
        ylim = axes[idx].get_ylim()
        low = min(xlim[0], ylim[0])
        high = max(xlim[1], ylim[1])
        line = np.linspace(low, high, 100)
        axes[idx].plot(line, line, 'k--', label='y = x')
        axes[idx].set_ylabel(f"true {xy_labels[idx]}")
        axes[idx].set_xlabel(f"converted {xy_labels[idx]}")
        axes[idx].set_title(f"{xy_labels[idx]}")
        if idx == 0:  # legend nur im ersten Plot
            axes[idx].legend(loc='upper left')
        axes[idx].grid(True)

        # wenn man mag kann man für die Übersichtlichkeit den Bildausschnitt noch anpassen
        axes[idx].set_xlim(right=ylim[1])

    # D/gamma extra
    # saving the RRMSE of D/gamma
    rrmse_ratio = []
    rsq_ratio = []
    for trial_set in range(len(labels)):
        if len(labels) == 1:  # sampling-Situation, in der D/gamma retrieved wurden ohne erst D und dann gamma zu holen
            temp = ordered_temp
            temp2 = D_over_gamma_converted

        else:  # situation in der man die Kohorten mit gleichem Alter kennt
            temp = [row[0] / row[1] for row in D_gamma_true_all_trials[trial_set]]
            temp2 = [row[0] / row[1] for row in D_gamma_converted_all_trials[trial_set]]

        axes[2].scatter(temp2, temp, marker=markers[trial_set], alpha=0.5, label=labels[trial_set])

        # compute the error
        squared_errors = [(temp[i] - temp2[i]) ** 2 for i in range(len(temp))]
        rmse = np.sqrt(np.mean(squared_errors))
        rrmse = compute_RRMSE(rmse, np.mean(temp))
        rsq = compute_R_square(np.sum(squared_errors), np.sum([(temp[i] - np.mean(temp)) ** 2 for i in range(len(temp))]))

        print(f'D/gamma:')
        print(f'Root mean squared error: {rmse}')
        print(f'Relative root mean squared error: {rrmse}')
        print(f'R² = {rsq}\n')
        rrmse_ratio.append(rrmse)
        rsq_ratio.append(rsq)

    # RRMSE
    mean_rrmse = np.mean(rrmse_ratio)
    std_rrmse = np.std(rrmse_ratio)

    # R²
    mean_rsq = np.mean(rsq_ratio)
    std_rsq = np.std(rsq_ratio)

    if len(labels) == 1:  # mean & std machen keinen Sinn bei len(labels) == 1
        axes[2].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f}', transform=axes[2].transAxes)
        axes[2].text(0.2, 0.05, f'R² = {mean_rsq:.3f}', transform=axes[2].transAxes)
    else:
        axes[2].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f} +/- {std_rrmse:.3f}', transform=axes[2].transAxes)
        axes[2].text(0.2, 0.05, f'R² = {mean_rsq:.3f} +/- {std_rsq:.3f}', transform=axes[2].transAxes)

    # draw the bisection line (first, identify domain & range, then plot)
    axes[2].set_xlim(left=0)
    axes[2].set_ylim(bottom=0)
    xlim = axes[2].get_xlim()
    ylim = axes[2].get_ylim()
    low = min(xlim[0], ylim[0])
    high = max(xlim[1], ylim[1])
    line = np.linspace(low, high, 100)
    axes[2].plot(line, line, 'k--', label='y = x')

    axes[2].set_ylabel("true D/gamma")
    axes[2].set_xlabel("converted D/gamma")
    axes[2].set_title("D/gamma")
    # axes[2].legend(loc='upper left')
    axes[2].grid(True)

    # wenn man mag kann man für die Übersichtlichkeit den Bildausschnitt noch anpassen
    # variant 1: gut für Fits mit kleinem Fehler
    # right_lim = max(ylim[1], xlim[1])
    # axes[2].set_xlim(right=right_lim)
    # axes[2].set_ylim(top=right_lim)
    # variant 2: gut für fits mit großem Fehler
    axes[2].set_xlim(right=ylim[1])

    plt.show()

    # D/gamma Plots: curve fit und parity (das sind die gleichen plots wie eben nur nochmal auf ner gemeinsamen Figure)
    if len(labels) == 1:
        fig, axes = plt.subplots(1, 2)
        fig.suptitle(rf't-cohorts unknown: Curve Fit and Parity Plot for $D/\gamma$')  #, R² > {thr_rsquared}')

        # curve fit plot
        axes[0].scatter(estimated_ratios, true_ratios, marker=markers[0], alpha=0.5, label=labels[0], c=colors)
        axes[0].plot(ordered_temp2, D_over_gamma_converted)

        cbar = plt.colorbar(sm, ax=axes[0])
        cbar.set_label('true D')  # Label for colorbar
        axes[0].set_ylabel("true D/gamma")
        axes[0].set_xlabel("estimated D/gamma")
        axes[0].set_title("Curve Fit")
        axes[0].legend()
        axes[0].grid(True)

        # parity plot
        axes[1].scatter(temp2, temp, marker=markers[0], alpha=0.5, label=labels[0])

        axes[1].text(0.2, 0.1, f'RRMSE = {mean_rrmse:.3f}', transform=axes[1].transAxes)
        axes[1].text(0.2, 0.05, f'R² = {mean_rsq:.3f}', transform=axes[1].transAxes)

        xlim = axes[1].get_xlim()
        ylim = axes[1].get_ylim()
        low = min(xlim[0], ylim[0])
        high = max(xlim[1], ylim[1])
        line = np.linspace(low, high, 100)
        axes[1].plot(line, line, 'k--', label='y = x')
        axes[1].set_ylabel("true D/gamma")
        axes[1].set_xlabel("converted D/gamma")
        axes[1].set_title("Parity Plot")
        axes[1].legend(loc='upper left')
        axes[1].grid(True)

        plt.show()

    # Fehler und Parameterwerte über t -- funktioniert nur wenn t-Kohorten bekannt
    if len(labels) > 1:
        # Fehler in D & gamma über t
        fig, axes = plt.subplots(2, 3)  # , figsize=(24, 8))
        fig.suptitle(f'Quality of fit statistics: fitting results vs. true parameter values over t')

        # wenn man das nicht vorher definiert hat
        # rrmse_ratio = []

        rrmses = [rrmse_D_array, rrmse_gamma_reconversion_array, rrmse_ratio, rrmse_gamma_array]
        rrmse_names = ['RRMSE of D', 'RRMSE of gamma 2nd iteration', 'RRMSE of D/gamma', 'RRMSE of gamma']

        rsqs = [rsq_D_array, rsq_gamma_reconversion_array, rsq_ratio, rsq_gamma_array]
        rsq_names = ['R² of D', 'R² of gamma 2nd iteration', 'R² of D/gamma', 'R² of gamma']
        for stat_idx, stat in enumerate([[rrmses, rrmse_names], [rsqs, rsq_names]]):
            for idx in [0, 1, 2]:
                axes[stat_idx, idx].plot(fixed_ts, stat[0][idx])
                mean = np.mean(stat[0][idx])
                std = np.std(stat[0][idx])
                axes[stat_idx, idx].text(0.05, 0.05, f'Avg. +/- 1 std.: {mean:.5f} +/- {std:.5f}', transform=axes[stat_idx, idx].transAxes)
                axes[stat_idx, idx].ticklabel_format(useOffset=False, style='plain')
                axes[stat_idx, idx].set_title(f'{stat[1][idx]}')
                axes[stat_idx, idx].set_ylabel(f'statistic')
                axes[stat_idx, idx].set_xlabel(f't')
                axes[stat_idx, idx].grid(True)
        plt.show()

        # Beziehung t <-> D
        # Plot der Parameterwerte a, b, c über t
        fig = plt.figure(figsize=(10, 6))
        gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1])  # 2 rows, 3 columns
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])
        axes = [ax1, ax2, ax3]
        ax4 = fig.add_subplot(gs[1, :])
        fig.suptitle(f'Parameters a, b, c from f(D^) = a sqrt(D^) + bD^ + c over t')
        params = ['a', 'b', 'c']

        for idx in range(3):
            axes[idx].plot(fixed_ts, [row[idx] for row in popts_D_array], marker='o', label=params[idx])
            axes[idx].grid(True)
            # axes[idx].ticklabel_format(useOffset=False, style='plain')
            axes[idx].set_ylabel(f'{params[idx]}')
            axes[idx].set_xlabel(f't')

            dparam_dt = np.gradient([row[idx] for row in popts_D_array], fixed_ts)
            ax4.plot(fixed_ts, dparam_dt, marker='x', label=f'd{params[idx]}/dt')
            ax4.set_ylabel('d/dt')
            ax4.set_xlabel(f't')
            ax4.legend()
            ax4.grid(True)
        plt.show()

        # Beziehung t <-> gamma
        if gamma_fit == 'pwlinear':
            # ich hab für jedes D ein eigenes a, soll ich da einfach den Durchschnitt nehmen?
            fig, axes = plt.subplots(1, 3)
            fig.suptitle(f'Parameters a, b, c from f(gamma^, D^) = a gamma^ + b f(D^) + c')

            for idx in range(3):
                param_vals_idx = [row[idx] for row in popts_gamma_array]
                param_means_to_plot = [np.mean(param_t) for param_t in param_vals_idx]
                axes[idx].plot(fixed_ts, param_means_to_plot, marker='o', label=f'means of {params[idx]}')
                axes[idx].set_ylabel(f'means of {params[idx]}')
                axes[idx].set_xlabel(f't')
                axes[idx].grid(True)

        else:
            fig = plt.figure()
            gs = gridspec.GridSpec(4, 3, height_ratios=[1, 1, 1, 1])  # 4 rows, 3 columns
            axes = []
            for row in range(3):
                ax1 = fig.add_subplot(gs[row, 0])
                ax2 = fig.add_subplot(gs[row, 1])
                ax3 = fig.add_subplot(gs[row, 2])
                axes.append([ax1, ax2, ax3])
            ax4 = fig.add_subplot(gs[3, :])

            fig.suptitle(f'Parameters m, n, o from i(D^) = m f(D^)² + n f(D^) + o, i in a, b, c')
            params_mno = ['m', 'n', 'o']
            colors = plt.cm.viridis.colors  # colormap damit fct & gradient dieselbe Farbe haben können

            gradients = []

            iteration = 0

            for idx1 in range(3):  # iterate through a, b, c
                for idx2 in range(3):  # iterate through m, n, o
                    color = colors[round(len(colors) * iteration / 9)]  # get the color
                    param_vals_idx = [t_mtx[idx1, idx2] for t_mtx in popts_gamma_quadr_array]
                    gradient = np.gradient(param_vals_idx, fixed_ts)
                    # axes[idx1][idx2].ticklabel_format(useOffset=False, style='plain')
                    axes[idx1][idx2].plot(fixed_ts, param_vals_idx, marker='o', color=color, alpha=0.7)
                    ax4.plot(fixed_ts, gradient, color=color, alpha=0.7,
                             label=f'd/dt {params_mno[idx2]} of {params[idx1]}')
                    axes[idx1][idx2].set_ylabel(f'{params_mno[idx2]} of {params[idx1]}')
                    axes[idx1][idx2].set_xlabel(f't')
                    axes[idx1][idx2].grid(True)

                    iteration += 1
            ax4.grid(True)
            ax4.legend()
        plt.show()

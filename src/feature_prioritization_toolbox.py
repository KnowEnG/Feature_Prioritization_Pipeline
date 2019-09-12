"""
@author: The KnowEnG dev team
"""
import os
import subprocess

import numpy as np
import pandas as pd

from sklearn.preprocessing import normalize

import knpackage.toolbox as kn
import knpackage.distributed_computing_utils as dstutil
import knpackage.data_cleanup_toolbox as datacln

EPSILON_0 = 1e-7

def run_correlation(run_parameters):
    """ perform feature prioritization

    Args:
        run_parameters: parameter set dictionary.
    """
    max_cpu = run_parameters["max_cpu"]
    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')

    phenotype_df = kn.get_spreadsheet_df(run_parameters["phenotype_name_full_path"])
    spreadsheet_df = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    phenotype_df = phenotype_df.T

    len_phenotype = len(phenotype_df.index)
    array_of_jobs = range(0, len_phenotype)

    for i in range(0, len_phenotype, max_cpu):
        jobs_id = array_of_jobs[i:i + max_cpu]
        number_of_jobs = len(jobs_id)

        zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, jobs_id)
        dstutil.parallelize_processes_locally(run_correlation_worker, zipped_arguments, number_of_jobs)
    write_phenotype_data_all(run_parameters)
    kn.remove_dir(run_parameters["results_tmp_directory"])

def run_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, job_id):
    """ core function for parallel run_correlation
    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        job_id:          parallel iteration number
    """
    # selects the ith row in phenotype_df

    np.random.seed(job_id)

    phenotype_df = phenotype_df.iloc[[job_id], :]

    spreadsheet_df, phenotype_df, msg = datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)

    pc_array = get_correlation(spreadsheet_df, phenotype_df, run_parameters)

    feature_name_list = spreadsheet_df.index
    phenotype_name = phenotype_df.index.values[0]
    generate_correlation_output(pc_array, phenotype_name, feature_name_list, run_parameters)

def generate_correlation_output(pc_array, phenotype_name, feature_name_list, run_parameters):
    """ Save final output of correlation
    
    Args:
        pc_array: pearson correlation coefficient array
        phenotype_name: name of the phenotype
        feature_name_list: list of the features correlated (size of pc_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """
    phenotype_name_list = np.repeat(phenotype_name, len(feature_name_list))
    baseline_score      = pc_array
    pc_array            = abs(pc_array)
    viz_score           = (pc_array - min(pc_array)) / (max(pc_array) - min(pc_array))
    pc_array            = np.round(pc_array,       8)
    viz_score           = np.round(viz_score,      8)
    baseline_score      = np.round(baseline_score, 8)
    
    output_val          = np.column_stack((phenotype_name_list, feature_name_list, pc_array, viz_score, baseline_score))
    df_header           = ['Response', 'Feature_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    
    result_df           = pd.DataFrame(columns=df_header)
    result_df['Response']                   = phenotype_name_list
    result_df['Feature_ID']                 = feature_name_list
    result_df['quantitative_sorting_score'] = pc_array
    result_df['visualization_score']        = viz_score
    result_df['baseline_score']             = baseline_score
    result_df = result_df.sort_values("visualization_score", ascending=0)
    result_df.index     = range(result_df.shape[0])

    write_one_phenotype(result_df, phenotype_name, feature_name_list, run_parameters)


def run_bootstrap_correlation(run_parameters):
    """ perform feature prioritization using bootstrap sampling

    Args:
        run_parameters: parameter set dictionary.
    """
    max_cpu             = run_parameters["max_cpu"]
    run_parameters["results_tmp_directory"] = kn.create_dir(run_parameters["results_directory"], 'tmp')

    phenotype_df        = kn.get_spreadsheet_df(run_parameters["phenotype_name_full_path"])
    spreadsheet_df      = kn.get_spreadsheet_df(run_parameters["spreadsheet_name_full_path"])
    phenotype_df        = phenotype_df.T
    n_bootstraps        = run_parameters["number_of_bootstraps"]

    len_phenotype = len(phenotype_df.index)
    array_of_jobs = range(0, len_phenotype)

    for i in range(0, len_phenotype, max_cpu):
        jobs_id = array_of_jobs[i:i + max_cpu]
        number_of_jobs = len(jobs_id)

        zipped_arguments = dstutil.zip_parameters(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, jobs_id)

        dstutil.parallelize_processes_locally(run_bootstrap_correlation_worker, zipped_arguments, number_of_jobs)
    write_phenotype_data_all(run_parameters)
    kn.remove_dir(run_parameters["results_tmp_directory"])


def run_bootstrap_correlation_worker(run_parameters, spreadsheet_df, phenotype_df, n_bootstraps, job_id):
    """  core function for parallel run_bootstrap_correlation

    Args:
        run_parameters:  dict of parameters
        spreadsheet_df:  spreadsheet data frame
        phenotype_df:    phenotype data frame
        n_bootstraps:    number of bootstrap samples to use
        job_id:          parallel iteration number
    """

    np.random.seed(job_id)

    phenotype_df      = phenotype_df.iloc[[job_id], :]

    spreadsheet_df, phenotype_df, msg = datacln.check_input_value_for_gene_prioritazion(spreadsheet_df, phenotype_df)

    pearson_array     = get_correlation(spreadsheet_df, phenotype_df, run_parameters)
    borda_count       = np.zeros(spreadsheet_df.shape[0])
    gm_accumulator    = np.ones(spreadsheet_df.shape[0])
    for bootstrap_number in range(0, n_bootstraps):
        # note: rows_sampling_fraction is not currently a run parameter
        spreadsheet_sample_df, phenotype_sample_df = sample_dfs(\
            spreadsheet_df, phenotype_df, 1.0, run_parameters["cols_sampling_fraction"])
        pc_array = get_correlation(\
            spreadsheet_sample_df, phenotype_sample_df, run_parameters)
        borda_count    = sum_array_ranking_to_borda_count(borda_count, np.abs(pc_array))
        gm_accumulator = (np.abs(pc_array) + EPSILON_0) * gm_accumulator

    pcc_gm_array       = gm_accumulator ** (1 / n_bootstraps)
    borda_count        = borda_count / n_bootstraps
    phenotype_name     = phenotype_df.index.values[0]
    feature_name_list  = spreadsheet_df.index
    viz_score          = (borda_count - min(borda_count)) / (max(borda_count) - min(borda_count))

    generate_bootstrap_correlation_output(borda_count, viz_score, pearson_array,
                                          phenotype_name, feature_name_list, run_parameters)


def generate_bootstrap_correlation_output(borda_count, viz_score, pearson_array, 
                                          phenotype_name, feature_name_list, run_parameters):
    """ Save final output of correlation

    Args:
        pearson_array: pearson correlation coefficient array
        phenotype_name: name of the phenotype
        feature_name_list: list of the features correlated (size of pearson_array
        run_parameters: dictionary of run parameters with key 'results_directory'
    """
    phenotype_name_list = np.repeat(phenotype_name, len(feature_name_list))
    viz_score           = np.round(viz_score, 8)
    borda_count         = np.round(borda_count, 8)
    pearson_array       = np.round(pearson_array, 8)

    output_val          = np.column_stack((phenotype_name_list, feature_name_list, borda_count, viz_score, pearson_array))
    df_header           = ['Response', 'Feature_ID', 'quantitative_sorting_score', 'visualization_score', 'baseline_score']
    
    result_df           = pd.DataFrame(columns=df_header)
    result_df['Response']                   = phenotype_name_list
    result_df['Feature_ID']                 = feature_name_list
    result_df['quantitative_sorting_score'] = borda_count
    result_df['visualization_score']        = viz_score
    result_df['baseline_score']             = pearson_array
    result_df = result_df.sort_values("visualization_score", ascending=0)
    result_df.index     = range(result_df.shape[0])

    write_one_phenotype(result_df, phenotype_name, feature_name_list, run_parameters)



def get_correlation(spreadsheet_df, phenotype_df, run_parameters):
    """ correlation function definition for all run methods

    Args:
        spreadsheet_df: features x samples
        phenotype_df: one x samples
        run_parameters: with key 'correlation_measure'

    Returns:
        correlation_array: features x one
    """

    if run_parameters['correlation_measure'] == 'pearson':

        spreadsheet_mat        = spreadsheet_df.values
        phenotype_arr          = phenotype_df.values[0]
        spreadsheet_mat        = spreadsheet_mat - spreadsheet_mat.mean(axis=1).reshape((-1, 1))
        phenotype_arr          = phenotype_arr - phenotype_arr.mean()
        spreadsheet_mat_var    = np.std(spreadsheet_mat, axis=1)
        phenotype_arr_var      = np.std(phenotype_arr)
        numerator              = spreadsheet_mat.dot(phenotype_arr)
        denominator            = spreadsheet_mat_var * phenotype_arr_var * spreadsheet_mat.shape[1]

        with np.errstate(divide='ignore', invalid='ignore'):
            correlation_array                 = np.true_divide(numerator, denominator)
            correlation_array[denominator==0] = 0

    elif run_parameters['correlation_measure'] == 't_test':

        spreadsheet_mat = spreadsheet_df.values
        phenotype_arr   = phenotype_df.values[0]
        a               = spreadsheet_mat[:, phenotype_arr!=0]
        b               = spreadsheet_mat[:, phenotype_arr==0]
        d               = np.mean(a, axis=1) - np.mean(b, axis=1)
        denom           = np.sqrt(np.var(a, axis=1, ddof=1)/a.shape[1] + np.var(b, axis=1, ddof=1)/b.shape[1])

        with np.errstate(divide='ignore', invalid='ignore'):
            correlation_array                  = np.divide(d, denom)
            correlation_array[np.isnan(denom)] = 0

    elif run_parameters['correlation_measure'] == 'edgeR':

        # attempted implementation with rpy2, but it would freeze with multiprocessing
        # next best thing: dump data to feather files and call Rscript

        filename_suffix = '_' + str(os.getpid()) + '.feather'
        spreadsheet_df_path = os.path.join(\
            run_parameters["results_tmp_directory"], 'spreadsheet_df' + filename_suffix)
        phenotype_df_path = os.path.join(\
            run_parameters["results_tmp_directory"], 'phenotype_df' + filename_suffix)
        output_df_path = os.path.join(\
            run_parameters["results_tmp_directory"], 'output_df' + filename_suffix))

        spreadsheet_df.to_feather(spreadsheet_df_path)
        phenotype_df.to_feather(phenotype_df_path)

        subprocess.run(["Rscript", "--vanilla", get_edgeR_script_path(), \
            spreadsheet_df_path, phenotype_df_path, output_df_path])

        correlations_matrix = pd.read_feather(output_df_path)

        # correlations_matrix will already be in the correct order
        correlation_array = correlations_matrix[:, 0]

    return correlation_array

def sum_array_ranking_to_borda_count(borda_count, corr_array):
    """ sum to borda count with a contigous array added to borda count
    Args:
        borda_count: the current borda count - same size as correlation array
        corr_array:  the correlation array to rank and add to the count
    Returns:
        borda_count: the ranking of the correlation array added to the input borda count
    """
    num_elem = borda_count.size

    # either assign (no duplicate case) or enumerate the correlation array
    if num_elem == (np.unique(corr_array)).size:
        borda_count[np.argsort(corr_array)] += np.int_(sorted(np.arange(0, corr_array.size) + 1))
        return borda_count

    # enumerate the borda vote
    borda_add     = np.zeros(num_elem)
    enum_value    = 1
    sort_order    = np.argsort(corr_array)
    current_value = corr_array[sort_order[0]]

    for k in range(0, num_elem):
        if corr_array[sort_order[k]] != current_value:
            enum_value   += 1
            current_value = corr_array[sort_order[k]]
        borda_add[sort_order[k]] = enum_value

    # scale to the number of elements in the array -- philosopical choice here --
    borda_add = borda_add + (num_elem - enum_value)

    return borda_count + borda_add

def sample_dfs(spreadsheet_df, phenotype_df, rows_fraction, cols_fraction):
    """Selects a random subset of the samples according to `cols_fraction` and a
    random subset of the features according to `rows_fraction`. Returns the
    portion of `spreadsheet_df` containing the selected samples and having all
    zeros for the unselected features. Returns the portion of `phenotype_df`
    containing the selected samples.

    Args:
        spreadsheet_df: feature x sample spreadsheet as a df.
        phenotype_df: one x sample phenotype data as a df.
        rows_fraction: the fraction of rows to keep in the returned sample - [0 : 1].
        cols_fraction: the fraction of cols to keep in the returned sample - [0 : 1].

    Returns:
        spreadsheet_sample_df: a df containing samples from `spreadsheet_df`.
        phenotype_sample_df: a df containing samples from `phenotype_df`.
    """
    # the following approach ensures results consistent with earlier versions
    # of the FP image (i.e., deliberately avoiding pandas shortcuts here)
    # for consistency of results, we also have to do the np.random.permutation
    # call for features prior to the call for samples

    # determine which rows are to be DISCARDED
    original_feature_count         = spreadsheet_df.shape[0]
    discarded_feature_count        = int(np.round(original_feature_count * (1 - rows_fraction)))
    discarded_features_permutation = np.random.permutation(original_feature_count)[0:discarded_feature_count]

    # determine which columns are to be KEPT
    original_sample_count          = spreadsheet_df.shape[1]
    kept_sample_count              = int(np.round(original_sample_count * cols_fraction))
    kept_sample_permutation        = np.random.permutation(original_sample_count)[0:kept_sample_count]

    # create new dfs with only the kept columns
    new_spreadsheet_df             = spreadsheet_df.iloc[:, kept_sample_permutation]
    new_phenotype_df               = phenotype_df.iloc[:, kept_sample_permutation]

    # fill discarded features with zeros
    if discarded_feature_count:
        # note: need to copy the spreadsheet to ensure we aren't modifying a
        # shared underlying representation of the data
        # TODO do downstream calculations really need discarded features set to
        # zero, or can we just exclude the discarded features' rows from the df?
        # if the latter, we could avoid the deep copy
        new_spreadsheet_df = new_spreadsheet_df.copy(deep=True)
        new_spreadsheet_df.iloc[discarded_features_permutation, :] = 0

    return new_spreadsheet_df, new_phenotype_df

def zscore_dataframe(features_by_sample_df):
    """ zscore by rows for features x samples dataframe

    Args:
        features_by_sample_df: zscore along rows for features x phenotypes dataframe

    Returns:
        spreadsheet_df: rows add up to zero, normalized to the mean and std deveiation
    """
    zscore_df = (features_by_sample_df.sub(features_by_sample_df.mean(axis=1), axis=0)).truediv(
                    np.maximum(features_by_sample_df.std(axis=1), 1e-12), axis=0)
    return zscore_df


def write_one_phenotype(result_df, phenotype_name, feature_name_list, run_parameters):
    """ write the phenotype output file to the results directory and the temporary directory files

    Args:
        result_df:
        phenotype_name:
        feature_name_list:
        run_parameters:
        
    Output:
        {phenotype}_{method}_{correlation_measure}_{timestamp}_viz.tsv
    """

    top_gamma_of_sort   = run_parameters['top_gamma_of_sort']

    result_df.to_csv(get_output_file_name(run_parameters, 'results_directory', phenotype_name, 'viz'), header=True, index=False, sep='\t', float_format="%g")

    download_result_df                 = pd.DataFrame(data=None, index=None, columns=[phenotype_name])
    download_result_df[phenotype_name] = result_df['Feature_ID']
    download_result_df.to_csv(
        get_output_file_name(run_parameters, 'results_tmp_directory', phenotype_name, 'download'), header=True, index=False, sep='\t', float_format="%g")

    top_features                       = download_result_df.values[: top_gamma_of_sort]
    update_orig_result_df              = pd.DataFrame(np.in1d(feature_name_list, top_features).astype(int), index=feature_name_list, columns=[phenotype_name])
    update_orig_result_df.to_csv(
        get_output_file_name(run_parameters, 'results_tmp_directory', phenotype_name, 'original'), header=True, index=True, sep='\t', float_format="%g")


def write_phenotype_data_all(run_parameters):
    """ Post Processing: writes rows as features, cols as phenotypes, data is feature in top n for the phenotype T or F.

    Args:
        run_parameters: with field: 'results_directory'
        
    Output:
        ranked_features_per_phenotype_{method}_{correlation_measure}_{timestamp}_download.tsv
        top_features_per_phenotype_{method}_{correlation_measure}_{timestamp}_download.tsv
    """  
    tmp_dir       = run_parameters["results_tmp_directory"]
    dirList       = sorted(os.listdir(tmp_dir))
    download_list = []
    original_list = []

    for fileName in dirList:
        if (fileName[-12:] == 'download.tsv'):
            download_list.append(fileName)
        if (fileName[-12:] == 'original.tsv'):
            original_list.append(fileName)

    if (len(download_list) == 0 or len(original_list) == 0):
        return

    StartFileName = os.path.join(tmp_dir, original_list[0])
    src_df        = pd.read_csv(StartFileName, sep='\t', header=0, index_col=0)
    index_list    = src_df.index.values

    all_phenotypes_download_df = pd.DataFrame(data=None, index=None)
    all_phenotypes_original_df = pd.DataFrame(data=None, index=index_list)


    for fileName in download_list:
        tFileName                                  = os.path.join(tmp_dir, fileName)
        src_df                                     = pd.read_csv(tFileName, sep='\t', header=0, index_col=None)
        phenotype_name                             = src_df.columns.values[0]
        all_phenotypes_download_df[phenotype_name] = src_df[phenotype_name]

    for fileName in original_list:
        tFileName                                  = os.path.join(tmp_dir, fileName)
        src_df                                     = pd.read_csv(tFileName, sep='\t', header=0, index_col=0)          
        phenotype_name                             = src_df.columns.values[0]
        all_phenotypes_original_df[phenotype_name] = src_df[phenotype_name]

    all_phenotypes_download_df.index = range(1, all_phenotypes_download_df.shape[0]+1)

    all_phenotypes_download_df.to_csv(
        get_output_file_name(run_parameters, 'results_directory', 'ranked_features_per_response', 'download'), header=True, index=True, sep='\t', float_format="%g")

    all_phenotypes_original_df.to_csv(
        get_output_file_name(run_parameters, 'results_directory', 'top_features_per_response', 'download'), header=True, index=True, sep='\t', float_format="%g")

def get_output_file_name(run_parameters, dir_name_key, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: dir_name_key, "method" and "correlation_measure"
        dir_name_key:   run_parameters dictionary key for the output directory
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before type_suffix
        type_suffix:    the file type extenstion (default 'tsv') without period character

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    output_file_name = os.path.join(run_parameters[dir_name_key], prefix_string + '_' +
                                    run_parameters['method'] + '_' + run_parameters["correlation_measure"])

    output_file_name = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix
    return output_file_name

def get_edgeR_script_path():
    """Returns the path to the R script that provides the edgeR implementation.

    Args:
        None.

    Returns:
        str: The path to the R script that provides the edgeR implementation.

    """
    # we expect to find the script in a parallel directory called "r_src"
    py_script_dir = os.path.dirname(os.path.realpath(__file__))
    parent_dir = os.path.dirname(py_script_dir)
    r_script_dir = os.path.join(parent_dir, "r_src")

    return os.path.join(r_script_dir, 'edgeR_fp.R')

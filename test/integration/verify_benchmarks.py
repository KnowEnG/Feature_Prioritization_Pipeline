"""
sobh@illinois.edu

"""

import os
import filecmp

verification_dir = '../data/verification'
results_dir      = '../test/run_dir/results'

def verify_benchmark(BENCHMARK_name,BENCHMARK_YML) :

    run_command  = 'python3 ../src/feature_prioritization.py -run_directory ./run_dir -run_file ' + BENCHMARK_YML
    os.system(run_command)

    All_files_in_results_dir = os.listdir(results_dir)

    for f in All_files_in_results_dir:
        if BENCHMARK_name in f :
            RESULT    = os.path.join(results_dir,      f             )
            BENCHMARK = os.path.join(verification_dir, BENCHMARK_name+'.tsv')
            if filecmp.cmp(RESULT, BENCHMARK) == True:
                print(BENCHMARK, 'PASS' )
            else:
                print(BENCHMARK, 'FAIL' )

def main():
    BENCHMARK = {'pearson'    : [ 
                                 'BENCHMARK_1_FP_pearson.yml'
                               , '17-AAG_correlation_pearson'
                               , 'ranked_features_per_phenotype_correlation_pearson'
                               , 'top_features_per_phenotype_correlation_pearson'
                               ] 
               ,'bootstrap_pearson'     : [  
                                 'BENCHMARK_2_FP_bootstrap_pearson.yml'
                               , '17-AAG_bootstrap_correlation_pearson'
                               , 'AEW541_bootstrap_correlation_pearson'
                               , 'AZD0530_bootstrap_correlation_pearson'
                               , 'AZD6244_bootstrap_correlation_pearson'
                               , 'Erlotinib_bootstrap_correlation_pearson'
                               , 'Irinotecan_bootstrap_correlation_pearson'
                               , 'L-685458_bootstrap_correlation_pearson'
                               , 'LBW242_bootstrap_correlation_pearson'
                               , 'Lapatinib_bootstrap_correlation_pearson'
                               , 'Nilotinib_bootstrap_correlation_pearson'
                               , 'Nutlin-3_bootstrap_correlation_pearson'
                               , 'PD-0325901_bootstrap_correlation_pearson'
                               , 'PD-0332991_bootstrap_correlation_pearson'
                               , 'PF2341066_bootstrap_correlation_pearson'
                               , 'PHA-665752_bootstrap_correlation_pearson'
                               , 'PLX4720_bootstrap_correlation_pearson'
                               , 'Paclitaxel_bootstrap_correlation_pearson'
                               , 'Panobinostat_bootstrap_correlation_pearson'
                               , 'RAF265_bootstrap_correlation_pearson'
                               , 'Sorafenib_bootstrap_correlation_pearson'
                               , 'TAE684_bootstrap_correlation_pearson'
                               , 'TKI258_bootstrap_correlation_pearson'
                               , 'Topotecan_bootstrap_correlation_pearson'
                               , 'ZD-6474_bootstrap_correlation_pearson'
                               , 'ranked_features_per_phenotype_bootstrap_correlation_pearson'
                               , 'top_features_per_phenotype_bootstrap_correlation_pearson'
                               ]
               ,'t_test': [  
                                 'BENCHMARK_3_FP_t_test.yml'
                               , '17-AAG_correlation_t_test'
                               , 'ranked_features_per_phenotype_correlation_t_test'
                               , 'top_features_per_phenotype_correlation_t_test'
                               ]
               ,'bootstra_t_test': [  
                                 'BENCHMARK_4_FP_bootstrap_t_test.yml'
                               , '17-AAG_bootstrap_correlation_t_test'
                               , 'AEW541_bootstrap_correlation_t_test'
                               , 'AZD0530_bootstrap_correlation_t_test'
                               , 'AZD6244_bootstrap_correlation_t_test'
                               , 'Erlotinib_bootstrap_correlation_t_test'
                               , 'Irinotecan_bootstrap_correlation_t_test'
                               , 'L-685458_bootstrap_correlation_t_test'
                               , 'LBW242_bootstrap_correlation_t_test'
                               , 'Lapatinib_bootstrap_correlation_t_test'
                               , 'Nilotinib_bootstrap_correlation_t_test'
                               , 'Nutlin-3_bootstrap_correlation_t_test'
                               , 'PD-0325901_bootstrap_correlation_t_test'
                               , 'PD-0332991_bootstrap_correlation_t_test'
                               , 'PF2341066_bootstrap_correlation_t_test'
                               , 'PHA-665752_bootstrap_correlation_t_test'
                               , 'PLX4720_bootstrap_correlation_t_test'
                               , 'Paclitaxel_bootstrap_correlation_t_test'
                               , 'Panobinostat_bootstrap_correlation_t_test'
                               , 'RAF265_bootstrap_correlation_t_test'
                               , 'Sorafenib_bootstrap_correlation_t_test'
                               , 'TAE684_bootstrap_correlation_t_test'
                               , 'TKI258_bootstrap_correlation_t_test'
                               , 'Topotecan_bootstrap_correlation_t_test'
                               , 'ZD-6474_bootstrap_correlation_t_test'
                               , 'ranked_features_per_phenotype_bootstrap_correlation_t_test'
                               , 'top_features_per_phenotype_bootstrap_correlation_t_test'
                               ]
                }

    os.system('make env_setup')
    for key in BENCHMARK.keys(): 
        BENCHMARK_list = BENCHMARK[key]
        BENCHMARK_YML  = BENCHMARK_list[0]
        for BENCHMARK_name in BENCHMARK_list[1:] :
            verify_benchmark(BENCHMARK_name,BENCHMARK_YML)
            os.system('rm ./run_dir/results/*')

if __name__ == "__main__":
    main()

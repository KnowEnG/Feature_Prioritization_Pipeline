"""
@author: The KnowEnG dev team
"""

def correlation(run_parameters):
    """ feature prioritization """
    from feature_prioritization_toolbox import run_correlation
    run_correlation(run_parameters)

def bootstrap_correlation(run_parameters):
    """ feature prioritization with bootstrap"""
    from feature_prioritization_toolbox import run_bootstrap_correlation
    run_bootstrap_correlation(run_parameters)

SELECT = {
    "correlation": correlation,
    "bootstrap_correlation": bootstrap_correlation }

def main():
    """
    This is the main function to perform feature prioritization.
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters

    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()

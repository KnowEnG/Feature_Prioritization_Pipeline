import unittest
from unittest import TestCase
import numpy as np
import pandas as pd

import feature_prioritization_toolbox as gpt

class TestGet_t_test_correlation(TestCase):
    """  predicted correlation vs feature_prioritization_toolbox.get_correlation for 'correlation_measure' = 't_test'  """
    def test_get_t_test_correlation(self):
        """ check t_test against hand calculated value, and check that return values are real numbers
        """
        n_test_rows = 25
        run_parameters = {'correlation_measure': 't_test'}
        predicted_correlation_member = 3
        predicted_t_value = 0.20751433916

        drug_response = np.array([1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0])
        spreadsheet_array = np.array(
            [0.1, 0.2, 0.3, 0.3, 0.2, 0.2, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.2, 0.3, 0.3, 0.2, 0.1, 0.1])

        n_cols = spreadsheet_array.size
        response_mat = np.zeros((n_test_rows,spreadsheet_array.size))
        spreadsheet_mat = np.zeros((n_test_rows,spreadsheet_array.size))

        predicted_correlation_array = []

        for k in range(0, n_test_rows):
            if k == predicted_correlation_member or np.random.random() > 0.5:
                P = np.random.permutation(spreadsheet_array.size)
                response_mat[k,:] = drug_response[P]
                spreadsheet_mat[k,:] = spreadsheet_array[P]
                predicted_correlation_array.append(k)
            else:
                response_mat[k,:] = drug_response
                spreadsheet_mat[k,:] = np.random.rand(n_cols)

        spreadsheet_mat[predicted_correlation_member,:] = spreadsheet_array
        spreadsheet_mat[-1,:] = 0

        spreadsheet_df = pd.DataFrame(data=spreadsheet_mat, \
            index=range(n_test_rows), columns=range(n_cols))
        response_df = pd.DataFrame(data=drug_response, \
            index=range(n_cols), columns=['response']).T
        corr_arr = gpt.get_correlation(spreadsheet_df, response_df, run_parameters)

        self.assertAlmostEqual(corr_arr[predicted_correlation_member], predicted_t_value, msg='t_test correlation error')
        for k in range(0, n_test_rows-2):
            self.assertTrue(np.isfinite(corr_arr[k]), msg='t_test Correlation Array is Not Finite')
        self.assertTrue(np.isnan(corr_arr[n_test_rows-1]))

if __name__ == '__main__':
    unittest.main()

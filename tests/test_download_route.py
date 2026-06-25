import unittest
from unittest.mock import patch
from app import app

class TestDownloadRoute(unittest.TestCase):
    def setUp(self):
        self.app = app.test_client()
        self.app.testing = True

    @patch('app.initialize_gee')
    @patch('utils.gee_handler.get_slope_from_dem')
    @patch('utils.gee_handler.get_landcover_at_point')
    @patch('utils.gee_handler.get_ndti_turbidity')
    def test_download_report_success(self, mock_ndti, mock_lc, mock_slope, mock_init):
        mock_init.return_value = True
        mock_slope.return_value = 0.005
        mock_lc.return_value = "Forest"
        mock_ndti.return_value = 0.1

        response = self.app.get('/download_report?lat=-12.0&lon=-77.0&d50=0.5')

        if response.status_code != 200:
            print(f"Response Data: {response.data}")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content_type, 'application/pdf')
        self.assertTrue(response.headers['Content-Disposition'].startswith('attachment; filename=informe_sedimentos_'))

if __name__ == '__main__':
    unittest.main()

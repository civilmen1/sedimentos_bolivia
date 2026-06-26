import unittest
import math
from unittest.mock import MagicMock, patch
import sys

# Mock ee before importing gee_handler
mock_ee = MagicMock()
sys.modules['ee'] = mock_ee

from utils.gee_handler import get_slope_from_dem, get_map_url

class TestGEEHandler(unittest.TestCase):
    @patch('utils.gee_handler.ee')
    def test_get_slope_mock(self, mock_ee):
        # Setup mock return values
        mock_img = MagicMock()
        mock_ee.Image.return_value = mock_img
        mock_ee.Terrain.slope.return_value = mock_img
        mock_img.reduceRegion.return_value.getInfo.return_value = {'slope': 10}

        # get_slope_from_dem returns slope as m/m = tan(angle), not radians
        slope = get_slope_from_dem(-12.0, -77.0)
        self.assertAlmostEqual(slope, math.tan(math.radians(10)), places=4)

    @patch('utils.gee_handler.ee')
    def test_get_map_url_mock(self, mock_ee):
        mock_img = MagicMock()
        mock_ee.Image.return_value = mock_img
        mock_ee.Terrain.slope.return_value = mock_img
        mock_img.getThumbURL.return_value = 'http://thumb/url'

        url = get_map_url(-12.0, -77.0, layer_type='slope')
        self.assertEqual(url, 'http://thumb/url')

if __name__ == '__main__':
    unittest.main()

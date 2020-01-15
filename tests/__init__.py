"""General vars for tests."""

from taurenmd.core import Path

datap = Path(
    Path(__file__).myparents(),
    'data',
    )

export_data_expected = Path(datap, 'export_data_expected.csv')
export_data_expected_2 = Path(datap, 'export_data_expected_2.csv')
toptest = Path(datap, 'pcnaA_frame0.pdb')
toptest_cif = Path(datap, 'pcnaA_frame0.cif')
trajtest = Path(datap, 'pcnaA.dcd')

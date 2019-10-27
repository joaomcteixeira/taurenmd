
import taurenmd.cli_imagemol as cli_imagemol  # noqa: F401
import taurenmd.cli_noSol as cli_noSol  # noqa: F401
import taurenmd.cli_report as cli_report  # noqa: F401
import taurenmd.cli_rmsd as cli_rmsd  # noqa: F401
import taurenmd.cli_template as cli_template  # noqa: F401
import taurenmd.cli_trajedit as cli_trajedit  # noqa: F401


def test_main():
    cli_template.main('dummy_topo', 'dummy_traj')

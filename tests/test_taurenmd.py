
import taurenmd.cli_template as cli_template
import taurenmd.cli_imagemol as cli_imagemol
import taurenmd.cli_noSol as cli_noSol
import taurenmd.cli_report as cli_report
import taurenmd.cli_rmsd as cli_rmsd
import taurenmd.cli_trajedit as cli_trajedit


def test_main():
    cli_template.main('dummy_topo', 'dummy_traj')

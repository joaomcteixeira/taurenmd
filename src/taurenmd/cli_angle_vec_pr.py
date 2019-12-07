"""
Vector projection angle calculator.

Calculates the angle between two vector projections in a reference plane.

Protocol:
    
    - Given three selections (sel1-3), calculates their individual
        centre of geometry (cog1-3), # esto seria las 3 subunidades de PCNA
    - Defines a reference plane based on cog1-3 (ref123), # el plane que corta PCNA longitudinalmente (medio donut para ti medio para mi los dos con agujero)
    - calculates the COG of the three cog1-3 (refCOG) # esto es el centro geometrico de los centros geometricos de las subunidades de pcna
    - given another selection (selA), defines a vector that goes from
        refCOG to the cog of selA.  # este es el vector que va desde el centro geometrico de pcna hasta el centro geometrico de la chain A. (es nuestro puntero de reloj)
    - defines the projection of vctAB in plane ref123. # define el vetor que es la proyecion del vetor anterior en el plano de pcna (esto deberia ser lo mismo en nuestro caso, pero por si acaso no lo es, porque estoy implementando una generalizacion)
    - For each frame: # para cada frame
    - Calculates the COG of the cog1-3 for the i-frame (iCOG)  # calcula el plano no PCNA
    - calculates vector refCOG-iCOG
    - subtracts above vector to  i-frame positions  # esto como hablamos ayer traslada las posiciones para las del frame referencia, es lo mismo que meter los vectores a 0,0 pero lo hago en un paso antes.
    - Defines the vector vctAB for the i-frame. # vector del centro de PCNA al centro de la cadena A
    - calculates the projection of vctABi on the ref123 plane.
    - calculates the angle between projections of vctABi and vctAB.

"""
import argparse

from bioplottemplates.plots import param

from taurenmd import log
from taurenmd.libs import libcalc, libcli, libmda, libutil
from taurenmd.logger import S, T


ap = libcli.CustomParser(
    description=__doc__,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    # formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    'topology',
    help='Topology file.',
    type=str,
    )

ap.add_argument(
    'trajectory',
    help=(
        'Trajectory files. If given, multiple trajectories will be'
        'contactenated by order.'
        ),
    nargs='+',
    )

ap.add_argument(
    '-l',
    '--plane-selection',
    help=(
        'Three selections which center of geometry defines the '
        'the plane. For example: \'segid A\' \'segidB \' \'segid C\'.'
        ' The angle will be calculated between this plane along the trajectory '
        ' to the reference frame.'
        ),
    nargs=3
    )

ap.add_argument(
    '-f',
    '--frame',
    help='Calc distances to a frame instead of relative along traj.',
    default=0,
    type=int,
    )

ap.add_argument(
    '-s',
    '--start',
    help='Start frame for slicing.',
    default=None,
    type=int,
    )

ap.add_argument(
    '-e',
    '--stop',
    help='Stop frame for slicing: exclusive',
    default=None,
    type=int,
    )

ap.add_argument(
    '-p',
    '--step',
    help='Step value for slicing',
    default=None,
    type=int,
    )

ap.add_argument(
    '-v',
    '--plotvars',
    help=(
        'Plot variables. '
        'Example: -v xlabel=frames ylabel=RMSD color=red.'
        ),
    nargs='*',
    action=libcli.ParamsToDict,
    )


def load_args():
    """Load user arguments."""
    cmd = ap.parse_args()
    return cmd


def maincli():
    cmd = load_args()
    main(**vars(cmd))
    return


def main(
        topology,
        trajectory,
        plane_selection,
        selectionA,
        frame=0,
        start=None,
        stop=None,
        step=None,
        plotvars=None,
        **kwargs
        ):
    log.info(T('calculating angles'))

    u = libmda.mda_load_universe(topology, *list(trajectory))

    frame_slice = libutil.frame_slice(
        start=start,
        stop=stop,
        step=step,
        )
    log.info(S('for slice {}', frame_slice))

    log.info(T('calculating plane eq. for first frame'))
    u.trajectory[0]
    cog1 = u.select_atoms(plane_selection[0]).center_of_geometry()
    cog2 = u.select_atoms(plane_selection[1]).center_of_geometry()
    cog3 = u.select_atoms(plane_selection[2]).center_of_geometry()

    ref_plane_123 = libcalc.calc_plane_eq(
        cog1,
        cog2,
        cog3,
        )
    log.info(S('the equation is {}x + {}y + {}z = {}', ra, rb, rc, rd))
    
    refCOG = u.select_atoms('({}) or ({}) or ({})'.format(*plane_selection))

    vecCOG_A = np.linalg.norm([
        refCOG,
        u.select_atoms(selectionA).center_of_geometry(),
        ])


    log.info(T('Calculating angles'))
    angles = []
    for fi, ts in zip(
            range(len(u.trajectory))[frame_slice],
            u.trajectory[frame_slice],
            ):

        point1 = u.select_atoms(plane_selection[0]).center_of_geometry()
        point2 = u.select_atoms(plane_selection[1]).center_of_geometry()
        point3 = u.select_atoms(plane_selection[2]).center_of_geometry()
        
        a, b, c, d = libcalc.calc_plane_eq(point1, point2, point3)

        try:
            angles.append(libcalc.calc_angle(ra, rb, rc, a, b, c))
        except ValueError:
            # happens when ValueError: math domain error
            # this is division by zero, means math.acos(d) is 0
            log.info(S('ValueError returned angle 0.0 found'))
            angles.append(0.0)
    
    log.info(S('calculated a total of {} angles.', len(angles)))

    if plotvars is None:
        plotvars = dict()
    
    if 'labels' not in plotvars:
        plotvars['labels'] = 'plane: {}'.format(' and '.join(plane_selection))

    log.info(T('plot params:'))
    for k, v in plotvars.items():
        log.info(S('{} = {!r}', k, v))

    param.plot(
        list(range(len(u.trajectory))[frame_slice]),
        angles,
        **plotvars,
        )

    log.info(S('done'))
    return


if __name__ == '__main__':
    maincli()

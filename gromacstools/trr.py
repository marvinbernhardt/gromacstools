import numpy as np
from MDAnalysis.lib.formats.libmdaxdr import TRRFile

def read_series(trr_filename, n_frames, readx=True, readv=False, readf=False):
    """Reads a trr trajectory. Returns a dictionary with positions,
velocities, boxes, steps, times, lambdas; all as series (last
array dimension is the frame number)."""
    # read first frame
    with TRRFile(trr_filename, 'r') as trr_file:
        frame = next(trr_file)
        n_atoms = len(frame.x)

    # prepare output series
    boxes = np.zeros((3, 3, n_frames))
    times = np.zeros(n_frames)
    steps = np.zeros(n_frames)
    lambdas = np.zeros(n_frames)
    if readx:
        positions = np.zeros((n_atoms*3, n_frames))
    else:
        positions = None
    if readv:
        velocities = np.zeros((n_atoms*3, n_frames))
    else:
        velocities = None
    if readf:
        forces = np.zeros((n_atoms*3, n_frames))
    else:
        forces = None

    # read n_frames frames
    with TRRFile(trr_filename, 'r') as trr_file:
        for frame_nr, frame in enumerate(trr_file):
            boxes[:,:,frame_nr] = frame.box
            times[frame_nr] = frame.time
            steps[frame_nr] = frame.step
            lambdas[frame_nr] = frame.lmbda
            if frame.hasx:
                positions[:, frame_nr] = frame.x.flat
            if frame.hasv:
                velocities[:, frame_nr] = frame.v.flat
            if frame.hasf:
                forces[:, frame_nr] = frame.f.flat

            if frame_nr == n_frames - 1:
                break
        # check if more frames read than available
        else:
            raise ValueError(f"n_frames is {n_frames} which is higher than the number of frames in {trr_filename}, which is {frame_nr+1}")

        series_dict = {'positions': positions, 'velocities': velocities, 'forces': forces, 'boxes': boxes, 'steps': steps, 'times': times, 'lambdas': lambdas}

        return series_dict

def write_series(trr_filename, series_dict):
    """Writes a trr trajectory. Takes a series dictionary
similar to what read_series returns."""

    def reshape_if_not_none(array, frame_nr):
        if array is None:
            return None
        else:
            return array[:, frame_nr].reshape((-1, 3))

    # write frames
    with TRRFile(trr_filename, 'w') as trr_file:
        for frame_nr in range(len(series_dict['steps'])):

            positions = reshape_if_not_none(series_dict['positions'], frame_nr)
            velocities = reshape_if_not_none(series_dict['velocities'], frame_nr)
            forces = reshape_if_not_none(series_dict['forces'], frame_nr)
            box = series_dict['boxes'][:,:,frame_nr]
            step = series_dict['steps'][frame_nr]
            time = series_dict['times'][frame_nr]
            lmbda = series_dict['lambdas'][frame_nr]
            n_atoms = len(positions)
            trr_file.write(positions, velocities, forces, box, step, time, lmbda, n_atoms)

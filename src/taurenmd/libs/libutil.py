def frame_slice(frame_id):
    if frame_id in ('all', None):
        return slice(None, None, None)

    elif isinstance(frame_id, list):
        return slice(*[int(i) for i in frame_id])

    else:
        raise NotImplementedError(f'Do not know what to do with: {frame_id}')

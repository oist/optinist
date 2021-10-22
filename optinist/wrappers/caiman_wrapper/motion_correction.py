def caiman_mc(info, opts=None):
    file_path = info['path']
    import numpy as np
    from caiman.source_extraction.cnmf.params import CNMFParams
    from caiman.motion_correction import MotionCorrect
    from caiman import load, save_memmap, load_memmap
    info = {}

    if opts is None:
        opts = CNMFParams()
    else:
        opts.change_params(params_dict=opts)

    mc = MotionCorrect(
        file_path, dview=None, **opts.get_group('motion'))

    mc.motion_correct(save_movie=True)
    border_to_0 = 0 if mc.border_nan == 'copy' else mc.border_to_0

    # memory mapping
    fname_new = save_memmap(
        mc.mmap_file, base_name='memmap_', order='C',border_to_0=border_to_0)

    # now load the file
    Yr, dims, T = load_memmap(fname_new)
    images = np.array(np.reshape(
        Yr.T, [T] + list(dims), order='F'))

    info['images'] = images

    return info


if __name__ == '__main__':
    import os
    info = {}
    file_path = os.path.join(
        '/Users', 'shogoakiyama', 'caiman_data', 
        'example_movies', 'Sue_2x_3000_40_-46.tif')
    info['caiman_mc'] = caiman_mc(file_path)
    print(info)

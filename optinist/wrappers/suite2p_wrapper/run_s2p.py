import suite2p

def run_s2p(file_path, opts=None):
    fnames = [file_path]

    if opts is None:
        from suite2p import default_ops
        opts = default_ops()

    data_path = []
    tiff_list = []
    for fname in fnames:
        data_path.append('/'.join(fname.split('/')[:-1]))
        tiff_list.append(fname.split('/')[-1])

    db = {
        'data_path': data_path,
        'save_dir': './logs',
        'save_path0': './logs',
        'tiff_list': tiff_list
    }

    info = suite2p.run_s2p(ops=opts, db=db)

    return info

if __name__ == '__main__':
    import os
    file_path = os.path.join(
        '/Users', 'shogoakiyama', 'Desktop', 'optinist', 
        'optinist', 'data', 'Sue_2x_3000_40_-46.tif')
    info = run_s2p(file_path)
    print(info.keys())

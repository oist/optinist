import os
import sys
sys.path.append('../optinist')

import json
import numpy as np
import pandas as pd
from PIL import Image


def run_code(wrapper_dict, flowList):
    print(flowList)
    print(wrapper_dict)

    info = {}
    inputs = None
    for item in flowList:
        print('-'*30)
        print('running')
        print(item)
        print('-'*30)
        if item.type == 'data':
            info[item.label] = {'path': item.path}
        elif item.type == 'algo':
            info[item.label] = wrapper_dict[item.label](info[prev_label])

        prev_label = item.label

    # save data to each direcotry
    from collections import OrderedDict
    results = OrderedDict()
    for item in flowList:
        results[item.label] = {}
        output = info[item.label]
        print(output.keys())

        if 'images' in output.keys():
            save_dir = os.path.join('./files', item.label, 'images')
            if not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)

            if len(output['images'].shape) == 2:
                output['images'] = output['images'][np.newaxis, :, :]

            for i in range(len(output['images'])):
                img = Image.fromarray(np.uint8(output['images'][i]))
                img.save(os.path.join(save_dir, f'{str(i)}.png'))

            results[item.label]['image_dir'] = save_dir

        if 'fluo' in output.keys():
            save_file = os.path.join('./files', item.label, 'fluo.json')

            pd.DataFrame(output['fluo']).to_json(save_file, indent=4)

            results[item.label]['fluo_path'] = save_file

    return {'message': 'success', 'outputPaths': results}


if __name__ == '__main__':
    run_code()

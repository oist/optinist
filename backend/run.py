import os
import sys
sys.path.append('../optinist')

import numpy as np
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
    result_path = []
    for item in flowList:
        output = info[item.label]

        if 'images' in output.keys():
            save_dir = os.path.join('./files', item.label, 'images')
            result_path.append(save_dir)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)

            if len(output['images'].shape) == 2:
                output['images'] = output['images'][np.newaxis, :, :]

            for i in range(len(output['images'])):
                img = Image.fromarray(np.uint8(output['images'][i]))
                img.save(os.path.join(save_dir, f'{str(i)}.png'))

    import random
    dummy_data = [{ "x": i, "y": random.uniform(100,0) } for i in range(0,20)]
    return {'message': 'success', "data": dummy_data, 'path': result_path}


if __name__ == '__main__':
    run_code()

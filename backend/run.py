import os
import sys
sys.path.append('../optinist')


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
            inputs = item.path
        elif item.type == 'algo':
            inputs = wrapper_dict[item.label](inputs)

    # info = {}
    # file_path = os.path.join(
    #     '/Users', 'shogoakiyama', 'caiman_data', 
    #     'example_movies', 'Sue_2x_3000_40_-46.tif')

    # from wrappers.caiman_wrapper import caiman_mc, caiman_cnmf, plot_contours_nb
    # info['caiman_mc'] = caiman_mc(file_path)
    # info['caiman_cnmf'] = caiman_cnmf(info['caiman_mc']['images'])

    return {'message': 'success'}

if __name__ == '__main__':
    run_code()

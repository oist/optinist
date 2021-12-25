import os
import pandas as pd
import cv2
import numpy as np
import json
from PIL import Image, ImageSequence
from wrappers import wrapper_dict
from wrappers.data_wrapper import *
from collections import OrderedDict
from typing import List


def get_dict_leaf_value(root_dict: dict, path_list: List[str]):
    path = path_list.pop(0)
    if(len(path_list) > 0):
        return get_dict_leaf_value(root_dict[path],path_list)
    else:
        return root_dict[path]


def save_tiff_to_json(tiff_file_path, maxidx=10):
    print("save_tiff_to_json")
    folder_path = os.path.dirname(tiff_file_path)
    file_name, ext = os.path.splitext(os.path.basename(tiff_file_path))

    # Tiff画像を読み込む
    tiffs = []
    image = Image.open(tiff_file_path)

    for i, page in enumerate(ImageSequence.Iterator(image)):

        if i >= maxidx:
            break

        page = np.array(page)

        if page.shape[-1] >= 200 and page.shape[-2] >= 200:
            page = cv2.resize(save_data, (200, 200))

        tiffs.append(page)

    tiffs = np.array(tiffs)

    images = []
    for i, _img in enumerate(tiffs):
        images.append(_img.tolist())

    pd.DataFrame(images).to_json(
        os.path.join(folder_path, f'{file_name}_{str(maxidx)}.json'), indent=4, orient="values")


def save_csv_to_json(csv_file_path):
    folder_path = os.path.dirname(csv_file_path)
    file_name, ext = os.path.splitext(os.path.basename(csv_file_path))
    pd.read_csv(csv_file_path).to_json(
        os.path.join(folder_path, f'{file_name}.json'), indent=4,orient="split")


def get_nest2dict(value):
    nwb_dict = {}

    for _k, _v in value.items():
        if _v['type'] == 'child':
            nwb_dict[_k] = _v['value']
        elif _v['type'] == 'parent':
            nwb_dict[_k] = get_nest2dict(_v['children'])

    return nwb_dict


def get_dict2nest(value):
    nwb_dict = {}
    for _k, _v in value.items():
        nwb_dict[_k] = {}
        if type(_v) is dict:
            nwb_dict[_k]['children'] = get_dict2nest(_v)
        else:
            nwb_dict[_k] = _v

    return nwb_dict


def string_to_float(params):
    for k, v in params.items():
        if isinstance(v, str) and v.isdecimal():
            v = float(v)
            if v.is_integer():
                v = int(v)
            params[k] = v

    return params


def get_algo_network(flowList):
    flowList = json.loads(flowList)
    nodeDict = {}
    startNodeList = []

    edgeList = flowList['elementListForRun']['edgeList']
    graph = {}

    # nodeを初期化
    for node in flowList['elementListForRun']['nodeList']:
        nodeDict[node['id']] = node
        graph[node['id']] = []
        if node['type'] != 'AlgorithmNode':
            startNodeList.append(node['id'])

    # 隣接リストを登録
    for edge in edgeList:
        graph[edge['source']].append(edge['target'])

    return graph, startNodeList, nodeDict


def run_algorithm(prev_info, item):
    if item['type'] == 'ImageFileNode':
        info = {'images': ImageData(item['data']['path'], '')}
    elif item['type'] == 'AlgorithmNode':
        # parameterをint, floatに変換
        if 'param' in item['data'].keys():
            params = string_to_float(item['data']['param'])
        else:
            params = {}

        wrapper = get_dict_leaf_value(
            wrapper_dict, item['data']['path'].split('/'))
        info = wrapper["function"](
            *prev_info.values(), params=params)
    else:
        assert False, 'run_algorithm error'

    return info


def get_results(info, item):
    results = OrderedDict()
    path = item['data']['path']
    results[path] = {}

    for k, v in info.items():
        if isinstance(v, ImageData):
            print("ImageData")
            results[path][k] = {}
            results[path][k]['path'] = v.json_path
            results[path][k]['type'] = 'images'
            results[path][k]['max_index'] = len(v.data)
        elif isinstance(v, TimeSeriesData):
            print("TimeSeriesData")
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'timeseries'
        elif isinstance(v, CorrelationData):
            print("CorrelationData")
            results[path][k] = {}
            results[path][k]['path'] = v.path
            results[path][k]['type'] = 'heatmap'
        else:
            pass

    return results

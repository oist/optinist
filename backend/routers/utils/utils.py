import gc
import numpy as np
import json
import copy


def nest2dict(value):
    nwb_dict = {}

    for _k, _v in value.items():
        if _v['type'] == 'child':
            nwb_dict[_k] = _v['value']
        elif _v['type'] == 'parent':
            nwb_dict[_k] = nest2dict(_v['children'])

    return nwb_dict


def dict2nest(value):
    nwb_dict = {}
    for _k, _v in value.items():
        nwb_dict[_k] = {}
        if type(_v) is dict:
            nwb_dict[_k]['children'] = dict2nest(_v)
        else:
            nwb_dict[_k] = _v

    return nwb_dict


# def string2float(params):
#     for k, v in params.items():
#         if isinstance(v, str) and v.isdecimal():
#             v = float(v)
#             if v.is_integer():
#                 v = int(v)
#             params[k] = v

#     return params


def check_types(params, default_params):
    for key in params.keys():
        if isinstance(params[key], dict):
            params[key] = check_types(params[key], default_params[key])
        else:
            if type(params[key]) != type(default_params[key]):
                data_type = type(default_params[key])
                p = params[key]
                if data_type == str:
                    params[key] = str(p)
                elif data_type == float:
                    params[key] = float(p)
                elif data_type == int:
                    params[key] = int(p)

    return params


def algo_network(flowList):
    flowList = json.loads(flowList)
    nodeDict = {}
    startNodeList = []

    edgeList = flowList['elementListForRun']['edgeList']
    graph = {}

    # nodeを初期化
    for node in flowList['elementListForRun']['nodeList']:
        nodeDict[node['id']] = node
        graph[node['id']] = []

    for node in flowList['elementListForRun']['nodeList']:
        if node['type'] == 'ImageFileNode':
            startNodeList.append(node['id'])

    for node in flowList['elementListForRun']['nodeList']:
        if node['type'] == 'CsvFileNode':
            startNodeList.append(node['id'])

    # 隣接リストを登録
    for edge in edgeList:
        graph[edge['source']] = {edge['target']: edge["targetHandle"].split("--")[1]}

    return graph, startNodeList, nodeDict

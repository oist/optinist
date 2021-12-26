from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import sys
import copy
import gc

from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException
from .utils import get_algo_network, run_algorithm, get_results

# NWB
from pynwb import NWBFile, NWBHDF5IO, get_manager
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries,
)
from .utils import get_nest2dict
from .nwb import nwb_config
from datetime import datetime
from dateutil.tz import tzlocal
import json
from .const import BASE_DIR
import tracemalloc

sys.path.append('../../optinist')

router = APIRouter()
manager = get_manager()


import linecache
import os
import tracemalloc

def display_top(snapshot, key_type='lineno', limit=10):
    snapshot = snapshot.filter_traces((
        tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
        tracemalloc.Filter(False, "<unknown>"),
    ))
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print("#%s: %s:%s: %.1f KiB" % (index, frame.filename, frame.lineno, stat.size / 1024))
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print('    %s' % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}


@router.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    # import run_pipeline

    await websocket.accept()
    # Wait for any message from the client

    message = await websocket.receive_text()
    graph, startNodeList, nodeDict = get_algo_network(message)

    message = json.loads(message)

    if message['nwbParam'] == {}:
        nwb_dict = copy.deepcopy(nwb_config)
    else:
        nwb_dict = get_nest2dict(message['nwbParam'])

    try:
        # import pdb; pdb.set_trace()
        tracemalloc.start()
        item = nodeDict[startNodeList[0]]
        prev_info = None
        info = None
        results = None

        while True:

            await websocket.send_json({
                'message': item['data']['label'] + ' started',
                'status': 'ready'
            })

            # run algorithm
            info = run_algorithm(prev_info, item)

            # NWB登録
            if item['type'] == 'ImageFileNode':
                nwb_dict['image_series']['external_file'] = info['images'].data
                nwbfile = nwb_add_acquisition(nwb_dict)
                nwbfile.create_processing_module(
                    name='ophys',
                    description='optical physiology processed data'
                )
                nwb_add_ophys(nwbfile)

                info['nwbfile'] = nwbfile

            # prev_infoを登録
            if 'nwbfile' not in info.keys():
                info['nwbfile'] = prev_info['nwbfile']
            else:
                # pass
                assert info is not None and 'nwbfile' in info.keys()

                # NWBを保存
                save_path = os.path.join(
                    BASE_DIR, item['data']['path'].split('.')[0].split('/')[-1])

                if not os.path.exists(save_path):
                    os.makedirs(save_path)

                # with NWBHDF5IO(
                #         os.path.join(save_path, 'tmp.nwb'), 'w',
                #         manager=manager) as f:
                #     f.write(info['nwbfile'])

            results = get_results(info, item)

            # Send message to the client
            await websocket.send_json({
                'message': item['data']['label'] + ' success',
                'status': 'success',
                'outputPaths': results
            })

            if len(graph[item['id']]) > 0:
                node_id = graph[item['id']][0]
                item = nodeDict[node_id]
            else:
                break

            prev_info = info

        with NWBHDF5IO(os.path.join(save_path, 'tmp.nwb'), 'w') as f:
            f.write(info['nwbfile'])

        del prev_info, info, nwbfile, results
        gc.collect()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics("lineno")

        print("---------------------------------------------------------")
        print("[ Top 10 ]")
        for stat in top_stats[:3]:
            print(stat)
        
        display_top(snapshot)

        await websocket.send_json({
            'message': 'completed',
            'status': 'completed'
        })

    except AlgorithmException as e:
        await websocket.send_json({
            'message': e.get_message(),
            'name': item['data']['label'],
            'status': 'error'
        })
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        await websocket.send_json({
            'message': message,
            'name': item['data']['label'],
            'status': 'error'
        })
    finally:
        gc.collect()
        print('Bye..')
        await websocket.close()

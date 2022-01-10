from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import os
import sys
import copy
import gc
import json
import tracemalloc

from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

# NWB
from pynwb import NWBFile, get_manager
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries,
)
from .utils import (
    get_nest2dict,
    save_nwb,
    get_algo_network,
    run_algorithm,
    get_results,
    display_top
)
from .nwb import nwb_config
from .const import BASE_DIR

sys.path.append('../../optinist')

router = APIRouter()
manager = get_manager()


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

            elif item['type'] == 'CsvFileNode':
                nwbfile = None
                info['nwbfile'] = nwbfile

            # prev_infoを登録
            if 'nwbfile' not in info.keys():
                info['nwbfile'] = prev_info['nwbfile']
            else:
                assert info is not None and 'nwbfile' in info.keys()

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

        # NWBを保存
        if info['nwbfile'] is not None:
            save_nwb(info['nwbfile'])

        del prev_info, info, nwbfile, results
        gc.collect()

        snapshot = tracemalloc.take_snapshot()
        top_stats = snapshot.statistics("lineno")
        display_top(snapshot, top_stats)

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

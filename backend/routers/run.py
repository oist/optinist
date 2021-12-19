from fastapi import APIRouter, WebSocket
from pydantic import BaseModel
from typing import Dict, List, Optional
import traceback
import sys
import json
from wrappers import wrapper_dict
from collections import OrderedDict
from wrappers.data_wrapper import *
from wrappers.optinist_exception import AlgorithmException

# NWB
from pynwb import NWBFile, NWBHDF5IO
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries,
)
from .utils import get_nest2dict
from .nwb import nwb_config
from datetime import datetime
from dateutil.tz import tzlocal

sys.path.append('../../optinist')

router = APIRouter()


class FlowItem(BaseModel):
    label: str
    path: Optional[str] = None
    type: str
    param: Optional[Dict] = {}


def get_dict_leaf_value(root_dict: dict, path_list: List[str]):
    path = path_list.pop(0)
    if(len(path_list) > 0):
        return get_dict_leaf_value(root_dict[path],path_list)
    else:
        return root_dict[path]


@router.get("/run/ready/{request_id}")
async  def run_ready(request_id: str):
    return {'message': 'ready...', 'status': 'ready', "requestId": request_id}


@router.websocket("/run")
async def websocket_endpoint(websocket: WebSocket):
    # import run_pipeline

    await websocket.accept()
    # Wait for any message from the client
    message = await websocket.receive_text()
    message = json.loads(message)
    flowList = list(map(lambda x: FlowItem(**x), message['nodeDataListForRun']))
    if message['nwbParam'] == {}:
        nwb_dict = nwb_config
    else:
        nwb_dict = get_nest2dict(message['nwbParam'])
    print('message: ', message)
    print('flowList: ', flowList)
    print('nwbParam: ', nwb_dict)

    try:
        for item in flowList:

            await websocket.send_json({'message': item.label+' started', 'status': 'ready'})

            # run algorithm
            info = None
            if item.type == 'image':
                info = {'path': ImageData(item.path, '')}
                # NWBを登録
                nwbfile = add_nwb_acquisition(nwb_dict)
            elif item.type == 'algo':
                wrapper = get_dict_leaf_value(wrapper_dict, item.path.split('/'))
                info = wrapper["function"](
                    *prev_info.values(), params=item.param)

            prev_info = info

            assert info is not None

            results = OrderedDict()
            results[item.path] = {}
            save_result(item, results, info)

            # NWBを保存
            print(item.path)
            save_path = item.path.split('.')[0]
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            with NWBHDF5IO(os.path.join(save_path, 'nwb.nwb'), 'w') as io:
                io.write(nwbfile)

            # Send message to the client
            await websocket.send_json({'message': item.label+' success', 'status': 'success', 'outputPaths': results})

        await websocket.send_json({'message': 'completed', 'status': 'completed'})
    except AlgorithmException as e:
        await websocket.send_json({'message': e.get_message(), 'name': item.label, 'status': 'error'})
    except Exception as e:
        message  = list(traceback.TracebackException.from_exception(e).format())[-1]
        print(traceback.format_exc())
        await websocket.send_json({'message': message, 'name': item.label, 'status': 'error'})
    finally:
        print('Bye..')
        await websocket.close()


def save_result(item, results, info):
    for k, v in info.items():
        if type(v) is ImageData:
            print("ImageData")
            results[item.path][k] = {}
            results[item.path][k]['path'] = v.path
            results[item.path][k]['type'] = 'images'
            results[item.path][k]['max_index'] = len(v.data)
        elif type(v) is TimeSeriesData:
            print("TimeSeriesData")
            results[item.path][k] = {}
            results[item.path][k]['path'] = v.path
            results[item.path][k]['type'] = 'timeseries'
        elif type(v) is CorrelationData:
            print("CorrelationData")
            results[item.path][k] = {}
            results[item.path][k]['path'] = v.path
            results[item.path][k]['type'] = 'heatmap'
        else:
            pass


def add_nwb_acquisition(nwb_dict):
    nwbfile = NWBFile(
        session_description=nwb_dict['session_description'],
        identifier=nwb_dict['identifier'],
        experiment_description=nwb_dict['experiment_description'],
        session_start_time=datetime.now(tzlocal()),
    )

    # 顕微鏡情報を登録
    device = nwbfile.create_device(
        name=nwb_dict['device']['name'], 
        description=nwb_dict['device']['description'],
        manufacturer=nwb_dict['device']['manufacturer']
    )

    # 光チャネルを登録
    optical_channel = OpticalChannel(
        name=nwb_dict['optical_channel']['name'], 
        description=nwb_dict['optical_channel']['description'], 
        emission_lambda=nwb_dict['optical_channel']['emission_lambda']
    )

    # imaging planeを追加
    nwbfile.create_imaging_plane(
        name=nwb_dict['imaging_plane']['name'],
        description=nwb_dict['imaging_plane']['description'],
        optical_channel=optical_channel,   # 光チャネル
        device=device,   # 電極デバイス
        imaging_rate=nwb_dict['imaging_plane']['imaging_rate'],   # 画像の比率Hz
        excitation_lambda=nwb_dict['imaging_plane']['excitation_lambda'], # 励起（れいき）波長
        indicator=nwb_dict['imaging_plane']['indicator'],   # カルシウムインディケーター
        location=nwb_dict['imaging_plane']['location'],
    )

    image_series = ImageSeries(
        name=nwb_dict['image_series']['name'],
        external_file=['image_file'],
        format='external',
        rate=nwb_dict['imaging_plane']['imaging_rate'],
        starting_frame=[nwb_dict['image_series']['starting_frame']]
    )

    nwbfile.add_acquisition(image_series)

    return nwbfile

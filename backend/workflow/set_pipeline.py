import yaml

from workflow.get_network import get_network
from workflow.params import get_snakemake_params, get_nwbfile
from workflow.set_file import set_imagefile, set_csvfile, set_algofile
from cui_api.snakemake import write_snakemake_config


def set_pipeline(runItem):
    rules_to_execute, last_outputs, all_outputs = get_workflow(runItem)

    flow_config = {
        "rules": rules_to_execute,
        "last_output": last_outputs,
    }
    
    write_snakemake_config(flow_config)


def get_workflow(runItem):
    # graph networkの解析
    nodeDict, edgeList, endNodeList = get_network(runItem)

    nwbfile = get_nwbfile(runItem.nwbParam)

    rules_to_execute = {}
    last_outputs = []
    all_outputs = {}

    for node in nodeDict.values():
        algo_label = node['data']['label']
        algo_path = node['data']['path']

        if node["type"] == 'ImageFileNode':
            set_imagefile(node, edgeList, nwbfile)
        elif node["type"] == "CsvFileNode":
            set_csvfile(node, edgeList)
        elif node["type"] == "AlgorithmNode":
            rule = set_algofile(node, edgeList, nodeDict)

            rules_to_execute[algo_label] = rule

            if node["id"] in endNodeList:
                last_outputs.append(rule["output"])

            all_outputs[rule["output"]] = {
                "label": algo_label,
                "path": algo_path,
            }

    return rules_to_execute, last_outputs, all_outputs


def dummy_run_pipeline(unique_id):
    import time
    i = 0
    os.makedirs(f"/tmp/optinist/{unique_id}")
    while True:
        print(f"i = {str(i)}")
        if i > 60:
            break

        if i == 5:
            with open(f"/tmp/optinist/{unique_id}/A.pkl", "wb") as f:
                info = {
                    'nodeId': 1,
                    'status': 'success',
                    'message': 'A success',
                    'name': 'A',
                    'outputPaths': {
                        "A_image": {
                            "path": f"/tmp/optinist/{unique_id}/A_images.json",
                            "type": "images",
                        },
                        "A_timeseries": {
                            "path": f"/tmp/optinist/{unique_id}/A_timeseries.json",
                            "type": "timeseries",
                        },
                    }
                }
                pickle.dump(info, f)
            

        if i == 10:
            with open(f"/tmp/optinist/{unique_id}/B.pkl", "wb") as f:
                info = {
                    'nodeId': 2,
                    'status': 'error',
                    'message': 'error reason',
                    'name': 'B',
                }
                pickle.dump(info, f)

        if i == 15:
            with open(f"/tmp/optinist/{unique_id}/C.pkl", "wb") as f:
                info = {
                    'nodeId': 3,
                    'status': 'success',
                    'message': 'C success',
                    'name': 'C',
                    'outputPaths': {
                        "C_image": {
                            "path": f"/tmp/optinist/{unique_id}/C_heatmp.json",
                            "type": "heatmap",
                        },
                        "C_timeseries": {
                            "path": f"/tmp/optinist/{unique_id}/C_timeseries.json",
                            "type": "timeseries",
                        },
                    }
                }
                pickle.dump(info, f)

        i += 1
        time.sleep(1)

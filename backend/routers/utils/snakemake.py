import os
import yaml


def create_snakemake_files(nodeDict):
    '''
    flowListを受け取り、Snakemakeに渡すconfig.yamlとして出力する。
    '''

    flow_config = {}
    rules_to_execute = {}
    prev_algo_output = None
    
    for i, item in nodeDict.items():
        if item['data']['type'] == 'input':
            algo_input = item["data"]["path"]
            algo_output = algo_input
        elif item["data"]["type"] == 'algorithm':
            algo_name = item["data"]["label"]
            algo_input = prev_output

            output_base_path = f"./files/{algo_name}"

            if not os.path.exists(output_base_path):
                print(f"Creating {output_base_path}")
                os.makedirs(output_base_path)

            algo_output = os.path.join(output_base_path, f"{algo_name}_out.pkl")

            rules_to_execute[algo_name] = {   
                "rule_file": f"rules/{algo_name}.smk",
                "input": algo_input,
                "param": item["data"]["param"],
                "output": algo_output
            }

        # 次のalgoのinputに渡すため、現在の出力ファイルのパスを保存
        prev_output = algo_output

    flow_config["rules"] = rules_to_execute
    flow_config["last_output"] = algo_output

    with open('./files/config.yaml', "w") as f:
        yaml.dump(flow_config, f)
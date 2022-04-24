# CUI execution
This document is how to execute on cui environment.

Set snakemake config yaml.
You can create config.yaml by frontend tool.
It created `config.yaml` in `/optinist/optinist` directory.
Previous `config.yaml` is located in `/tmp/optinist/${unique_id}/config.yaml`, so put into `/optinist/optinist` directory.


Set optinist environment.
```
conda activate optinist
```

move to optinist directory.
```
cd ${YOUR WORKING DIRECTORY}/optinist/optinist
```

run command
```
snakemake --cores 2 --forceall
```
import React from 'react'
import { useSelector } from 'react-redux'

import FormControl from '@mui/material/FormControl'
import MenuItem from '@mui/material/MenuItem'
import InputLabel from '@mui/material/InputLabel'
import FormHelperText from '@mui/material/FormHelperText'
import Select from '@mui/material/Select'
import ListSubheader from '@mui/material/ListSubheader'

import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from 'store/slice/DisplayData/DisplayDataType'
import { FILE_TYPE } from 'store/slice/InputNode/InputNodeType'
import { RootState } from 'store/store'
import { selectInputNode } from 'store/slice/InputNode/InputNodeSelectors'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'
import {
  selectPipelineLatestUid,
  selectPipelineNodeResultSuccessList,
} from 'store/slice/Pipeline/PipelineSelectors'
import { getFileName } from 'store/slice/FlowElement/FlowElementUtils'

export const FilePathSelect: React.FC<{
  dataType?: DATA_TYPE
  selectedNodeId: string | null
  selectedFilePath: string | null
  onSelect: (nodeId: string, filePath: string, dataType: DATA_TYPE, outputKey?: string) => void
  label?: string
}> = ({ dataType, selectedNodeId, selectedFilePath, onSelect, label }) => {
  const inputNodeFilePathInfoList = useSelector(
    (state: RootState) => {
      const inputNodes = selectInputNode(state)
      return Object.entries(inputNodes)
        .map(([nodeId, inputNode]) => ({
          nodeId,
          filePath: inputNode.selectedFilePath,
          fileType: inputNode.fileType,
          dataType: toDataTypeFromFileType(inputNode.fileType),
          nodeName: selectNodeLabelById(nodeId)(state),
        }))
        .filter(({ filePath }) => filePath != null)
        .filter(({ dataType: inputNodeDataType }) =>
          dataType != null ? inputNodeDataType === dataType : true,
        )
    },
    // todo 比較関数
  )

  const latestUid = useSelector(selectPipelineLatestUid)

  const algorithmNodeOutputPathInfoList = useSelector((state: RootState) => {
    if (latestUid != null) {
      const runResult = selectPipelineNodeResultSuccessList(state)
      return runResult.map(({ nodeId, nodeResult }) => {
        return {
          nodeId,
          nodeName: selectNodeLabelById(nodeId)(state),
          paths: Object.entries(nodeResult.outputPaths)
            .map(([outputKey, value]) => {
              return {
                outputKey,
                filePath: value.path,
                type: value.type,
              }
            })
            .filter(({ type }) =>
              dataType != null ? type === dataType : true,
            ),
        }
      })
    } else {
      return []
    }
  })

  const [open, setOpen] = React.useState(false)
  const handleClose = () => {
    setOpen(false)
  }

  const handleOpen = () => {
    setOpen(true)
  }

  const onSelectHandle = (
    nodeId: string,
    filePath: string,
    dataType: DATA_TYPE,
    outputKey?: string
  ) => {
    onSelect(nodeId, filePath, dataType, outputKey)
    handleClose()
  }

  const menuItemList: React.ReactElement[] = []
  inputNodeFilePathInfoList.forEach((pathInfo) => {
    const filePath = pathInfo.filePath
    if (Array.isArray(filePath)) {
      filePath.forEach((pathElm) => {
        menuItemList.push(
          <MenuItem
            value={`${pathInfo.nodeId}/${pathElm}`}
            onClick={() =>
              onSelectHandle(pathInfo.nodeId, pathElm ?? '', pathInfo.dataType)
            }
            key={pathInfo.nodeId}
          >
            {getFileName(pathElm)}
          </MenuItem>,
        )
      })
    } else {
      menuItemList.push(
        <MenuItem
          value={`${pathInfo.nodeId}/${pathInfo.filePath}`}
          onClick={() =>
            onSelectHandle(pathInfo.nodeId, filePath ?? '', pathInfo.dataType)
          }
          key={pathInfo.nodeId}
        >
          {pathInfo.nodeName}
        </MenuItem>,
      )
    }
  })
  algorithmNodeOutputPathInfoList.forEach((pathInfo) => {
    menuItemList.push(<ListSubheader>{pathInfo.nodeName}</ListSubheader>)
    pathInfo.paths.forEach((outputPath, i) => {
      menuItemList.push(
        <MenuItem
          value={`${pathInfo.nodeId}/${outputPath.filePath}`}
          onClick={() =>
            onSelectHandle(
              pathInfo.nodeId,
              outputPath.filePath,
              outputPath.type,
              outputPath.outputKey
            )
          }
          key={`${pathInfo.nodeId}/${outputPath.filePath}`}
        >
          {outputPath.outputKey}
        </MenuItem>,
      )
    })
  })

  return (
    <FormControl style={{ minWidth: 150, maxWidth: 220 }} variant="standard">
      <InputLabel>{!!label ? label : 'Select Item'}</InputLabel>
      <Select
        value={`${selectedNodeId}/${selectedFilePath}`}
        open={open}
        onClose={handleClose}
        onOpen={handleOpen}
      >
        {menuItemList}
      </Select>
      {inputNodeFilePathInfoList.length +
        algorithmNodeOutputPathInfoList.length ===
        0 && <FormHelperText error={true}>no data</FormHelperText>}
    </FormControl>
  )
}

function toDataTypeFromFileType(fileType: FILE_TYPE) {
  switch (fileType) {
    case FILE_TYPE_SET.IMAGE:
      return DATA_TYPE_SET.IMAGE
    case FILE_TYPE_SET.CSV:
      return DATA_TYPE_SET.CSV
    case FILE_TYPE_SET.HDF5:
      return DATA_TYPE_SET.HDF5
    case FILE_TYPE_SET.FLUO:
      return DATA_TYPE_SET.FLUO
    case FILE_TYPE_SET.BEHAVIOR:
      return DATA_TYPE_SET.BEHAVIOR
  }
}

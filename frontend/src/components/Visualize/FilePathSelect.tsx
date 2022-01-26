import React from 'react'
import { useSelector } from 'react-redux'

import FormControl from '@material-ui/core/FormControl'
import MenuItem from '@material-ui/core/MenuItem'
import InputLabel from '@material-ui/core/InputLabel'
import FormHelperText from '@material-ui/core/FormHelperText'
import Select from '@material-ui/core/Select'
import ListSubheader from '@material-ui/core/ListSubheader'

import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import { RootState } from 'store/store'
import { selectInputNode } from 'store/slice/InputNode/InputNodeSelectors'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import { selectAlgorithmNode } from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { selectOutputPaths } from 'store/slice/RunPipelineResult/RunPipelineResultSelectors'
import { toDataType } from 'store/slice/DisplayData/DisplayDataUtils'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'

export const FilePathSelect: React.FC<{
  dataType: string
  selectedNodeId: string | null
  selectedFilePath: string | null
  onSelect: (nodeId: string, filePath: string) => void
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
          nodeName: selectNodeLabelById(nodeId)(state),
        }))
        .filter(({ filePath }) => filePath != null)
        .filter(({ fileType }) => {
          switch (dataType) {
            case DATA_TYPE_SET.IMAGE:
              return fileType === FILE_TYPE_SET.IMAGE
            case DATA_TYPE_SET.CSV:
              return fileType === FILE_TYPE_SET.CSV
            default:
              return false
          }
        })
    },
    // todo 比較関数
  )

  const algorithmNodeOutputPathInfoList = useSelector((state: RootState) => {
    const algorithms = selectAlgorithmNode(state)
    const outputPaths = selectOutputPaths(state)
    if (outputPaths != null) {
      return Object.entries(algorithms)
        .filter(([nodeId, algoNode]) =>
          Object.keys(outputPaths).includes(algoNode.functionPath),
        )
        .map(([nodeId, algoNode]) => {
          const paths = Object.entries(outputPaths[algoNode.functionPath])
            .map(([outputKey, outputPath]) => ({
              outputKey,
              filePath: outputPath.path,
              type: toDataType(outputPath.type),
            }))
            .filter(({ type }) => type === dataType)
          return {
            nodeName: selectNodeLabelById(nodeId)(state),
            nodeId,
            paths,
          }
        })
        .filter(({ paths }) => paths.length > 0)
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

  const onSelectHandle = (nodeId: string, filePath: string) => {
    onSelect(nodeId, filePath)
    handleClose()
  }

  const menuItemList: React.ReactElement[] = []
  inputNodeFilePathInfoList.forEach((pathInfo) => {
    menuItemList.push(
      <MenuItem
        value={`${pathInfo.nodeId}/${pathInfo.filePath}`}
        onClick={() => onSelectHandle(pathInfo.nodeId, pathInfo.filePath ?? '')}
        key={pathInfo.nodeId}
      >
        {pathInfo.nodeName}
      </MenuItem>,
    )
  })
  algorithmNodeOutputPathInfoList.forEach((pathInfo) => {
    menuItemList.push(<ListSubheader>{pathInfo.nodeName}</ListSubheader>)
    pathInfo.paths.forEach((outputPath, i) => {
      menuItemList.push(
        <MenuItem
          value={`${pathInfo.nodeId}/${outputPath.filePath}`}
          onClick={() => onSelectHandle(pathInfo.nodeId, outputPath.filePath)}
          key={`${pathInfo.nodeId}/${outputPath.filePath}`}
        >
          {outputPath.outputKey}
        </MenuItem>,
      )
    })
  })
  return (
    <FormControl style={{ minWidth: 150, maxWidth: 220 }}>
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

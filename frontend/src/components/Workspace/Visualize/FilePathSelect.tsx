import { FC, useState, ReactElement } from "react"
import { useSelector } from "react-redux"

import { Divider } from "@mui/material"
import FormControl from "@mui/material/FormControl"
import FormHelperText from "@mui/material/FormHelperText"
import InputLabel from "@mui/material/InputLabel"
import ListSubheader from "@mui/material/ListSubheader"
import MenuItem from "@mui/material/MenuItem"
import Select from "@mui/material/Select"

import { FileNodeFactory } from "factories/FileNodeFactory"
import { DATA_TYPE } from "store/slice/DisplayData/DisplayDataType"
import { selectNodeLabelById } from "store/slice/FlowElement/FlowElementSelectors"
import { getFileName } from "store/slice/FlowElement/FlowElementUtils"
import { selectInputNode } from "store/slice/InputNode/InputNodeSelectors"
import {
  selectPipelineLatestUid,
  selectPipelineNodeResultSuccessList,
} from "store/slice/Pipeline/PipelineSelectors"
import { RootState } from "store/store"

export const FilePathSelect: FC<{
  dataType?: DATA_TYPE
  selectedNodeId: string | null
  selectedFilePath: string | null
  onSelect: (
    nodeId: string,
    filePath: string,
    dataType: DATA_TYPE,
    outputKey?: string,
  ) => void
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
          dataType: FileNodeFactory.toDataTypeFromFileType(inputNode.fileType),
          nodeName: selectNodeLabelById(nodeId)(state),
        }))
        .filter(({ filePath }) => filePath != null)
        .filter(({ dataType: inputNodeDataType }) =>
          dataType != null ? inputNodeDataType === dataType : true,
        )
    },
    (prevFilePathInfoList, nextFilePathInfoList) =>
      prevFilePathInfoList.length === nextFilePathInfoList.length &&
      prevFilePathInfoList.every(
        (prevInfo, index) =>
          prevInfo.filePath === nextFilePathInfoList[index].filePath &&
          prevInfo.fileType === nextFilePathInfoList[index].fileType &&
          prevInfo.dataType === nextFilePathInfoList[index].dataType &&
          prevInfo.nodeName === nextFilePathInfoList[index].nodeName,
      ),
  )

  const latestUid = useSelector(selectPipelineLatestUid)

  const algorithmNodeOutputPathInfoList = useSelector(
    (state: RootState) => {
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
    },
    (prevPathInfoList, nextPathInfoList) =>
      prevPathInfoList.length === nextPathInfoList.length &&
      prevPathInfoList.every((prevInfo, index) => {
        const nextInfo = nextPathInfoList[index]
        return (
          prevInfo.nodeId === nextInfo.nodeId &&
          prevInfo.nodeName === nextInfo.nodeName &&
          prevInfo.paths.length === nextInfo.paths.length &&
          prevInfo.paths.every(
            (prevPath, pathIndex) =>
              prevPath.outputKey === nextInfo.paths[pathIndex].outputKey &&
              prevPath.filePath === nextInfo.paths[pathIndex].filePath &&
              prevPath.type === nextInfo.paths[pathIndex].type,
          )
        )
      }),
  )

  const [open, setOpen] = useState(false)
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
    outputKey?: string,
  ) => {
    onSelect(nodeId, filePath, dataType, outputKey)
    handleClose()
  }

  const menuItemList: ReactElement[] = []
  inputNodeFilePathInfoList.forEach((pathInfo) => {
    const filePath = pathInfo.filePath
    if (Array.isArray(filePath)) {
      filePath.forEach((pathElm) => {
        menuItemList.push(
          <MenuItem
            value={
              pathInfo.nodeId && pathElm ? `${pathInfo.nodeId}/${pathElm}` : ""
            }
            onClick={() =>
              onSelectHandle(pathInfo.nodeId, pathElm ?? "", pathInfo.dataType)
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
          value={
            pathInfo.nodeId && pathInfo.filePath
              ? `${pathInfo.nodeId}/${pathInfo.filePath}`
              : ""
          }
          onClick={() =>
            onSelectHandle(pathInfo.nodeId, filePath ?? "", pathInfo.dataType)
          }
          key={pathInfo.nodeId}
        >
          {pathInfo.nodeName}
        </MenuItem>,
      )
    }
  })
  algorithmNodeOutputPathInfoList.forEach((pathInfo, index) => {
    menuItemList.push(
      <ListSubheader key={index}>
        <Divider textAlign="center">{pathInfo.nodeName}</Divider>
      </ListSubheader>,
    )
    pathInfo.paths.forEach((outputPath) => {
      menuItemList.push(
        <MenuItem
          value={
            pathInfo.nodeId && outputPath.filePath
              ? `${pathInfo.nodeId}/${outputPath.filePath}`
              : ""
          }
          onClick={() =>
            onSelectHandle(
              pathInfo.nodeId,
              outputPath.filePath,
              outputPath.type,
              outputPath.outputKey,
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
      <InputLabel>{label ? label : "Select Item"}</InputLabel>
      <Select
        value={
          selectedNodeId && selectedFilePath
            ? `${selectedNodeId}/${selectedFilePath}`
            : ""
        }
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

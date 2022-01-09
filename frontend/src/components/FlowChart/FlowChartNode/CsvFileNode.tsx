import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@material-ui/core/styles'
import { Typography } from '@material-ui/core'

import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectInputNodeDefined,
  selectInputNodeSelectedFilePath,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeSlice'
import { useDisplayDataTabAciton } from 'components/flextlayout/FlexLayoutHook'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'

const sourceHandleStyle: CSSProperties = {
  width: 8,
  height: 15,
  top: 15,
  border: '1px solid',
  borderColor: '#555',
  borderRadius: 0,
}

export const CsvFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <CsvFileNodeImple {...element} />
  } else {
    return null
  }
})

const CsvFileNodeImple = React.memo<NodeProps>(({ id: nodeId, selected }) => {
  const { setDisplayTab } = useDisplayDataTabAciton(nodeId)
  const dispatch = useDispatch()
  const filePath = useSelector(selectInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }
  const onClick = () => {
    if (filePath != null) {
      setDisplayTab(filePath, DATA_TYPE_SET.TABLE)
    }
  }
  const theme = useTheme()
  return (
    <div
      style={{
        height: '100%',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
      onClick={onClick}
    >
      <CsvFileSelect
        nodeId={nodeId}
        onChangeFilePath={onChangeFilePath}
        filePath={filePath ?? ''}
      />
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, 'csv', 'TableData')}
        style={sourceHandleStyle}
      />
    </div>
  )
})

const CsvFileSelect = React.memo<{
  nodeId: string
  filePath: string
  onChangeFilePath: (path: string) => void
}>(({ nodeId, filePath, onChangeFilePath }) => {
  const { setDisplayTab, deleteDisplayTab } = useDisplayDataTabAciton(nodeId)
  const {
    filePath: uploadedFilePath,
    onUploadFile,
    pending,
    uninitialized,
    progress,
    error,
  } = useFileUploader(FILE_TYPE_SET.CSV)
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
    if (uploadedFilePath != null) {
      onChangeFilePath(uploadedFilePath)
      setDisplayTab(uploadedFilePath, DATA_TYPE_SET.TABLE)
      if (filePath !== uploadedFilePath) {
        deleteDisplayTab(filePath) // 前回のfilePathを表示するtabは削除
      }
    }
  }
  const onSelectFile = (selectedPath: string) => {
    onChangeFilePath(selectedPath)
    setDisplayTab(selectedPath, DATA_TYPE_SET.TABLE)
    if (selectedPath !== filePath) {
      deleteDisplayTab(filePath) // 前回のfilePathを表示するtabは削除
    }
  }
  return (
    <>
      {!uninitialized && pending && progress != null && (
        <div style={{ marginLeft: 2, marginRight: 2 }}>
          <LinearProgressWithLabel value={progress} />
        </div>
      )}
      <FileSelect
        filePath={filePath}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFileHandle}
        fileTreeType={FILE_TREE_TYPE_SET.CSV}
        selectButtonLabel="Select CSV"
      />
      {error != null && (
        <Typography variant="caption" color="error">
          {error}
        </Typography>
      )}
    </>
  )
})

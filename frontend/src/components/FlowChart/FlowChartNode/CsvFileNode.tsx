import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@material-ui/core/styles'
import { Typography, IconButton } from '@material-ui/core'
import CloseOutlinedIcon from '@material-ui/icons/CloseOutlined'

import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectInputNodeDefined,
  selectInputNodeSelectedFilePath,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeSlice'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'
import {
  deleteFlowElementsById,
  edifFlowElementsLabelById,
} from 'store/slice/FlowElement/FlowElementSlice'

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
  const dispatch = useDispatch()
  const filePath = useSelector(selectInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
    const fileName = path.split('/').reverse()[0]
    dispatch(
      edifFlowElementsLabelById({
        nodeId,
        fileName,
      }),
    )
  }
  const theme = useTheme()

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowElementsById(nodeId))
  }

  return (
    <div
      style={{
        height: '100%',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
    >
      <IconButton
        aria-label="delete"
        style={{ color: 'black', position: 'absolute', top: -20, right: -5 }}
        onClick={onClickDeleteIcon}
      >
        <CloseOutlinedIcon />
      </IconButton>
      <CsvFileSelect
        nodeId={nodeId}
        onChangeFilePath={onChangeFilePath}
        filePath={filePath ? filePath.split('/').reverse()[0] : ''}
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
    }
  }
  const onSelectFile = (selectedPath: string) => {
    onChangeFilePath(selectedPath)
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

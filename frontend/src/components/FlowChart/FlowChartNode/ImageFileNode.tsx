import React, { CSSProperties } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@material-ui/core/styles'
import { Typography, IconButton } from '@material-ui/core'
import CloseOutlinedIcon from '@material-ui/icons/CloseOutlined'

import { FILE_TREE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import {
  selectInputNodeDefined,
  selectInputNodeSelectedFilePath,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setInputImageNodeFile } from 'store/slice/InputNode/InputNodeSlice'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'

import { useHandleColor } from './HandleColorHook'
import { FileSelect } from './FileSelect'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'
import { toHandleId, isValidConnection } from './FlowChartUtils'
import {
  deleteFlowElementsById,
  edifFlowElementsLabelById,
} from 'store/slice/FlowElement/FlowElementSlice'

// Connection部分のレイアウト
const sourceHandleStyle: CSSProperties = {
  width: '4%',
  height: '13%',
  top: 15,
  border: '1px solid',
  // borderColor: '#555',
  borderRadius: 0,
}

export const ImageFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <ImageFileNodeImple {...element} />
  } else {
    return null
  }
})

const ImageFileNodeImple = React.memo<NodeProps>(
  ({ id: nodeId, selected: elementSelected }) => {
    const dispatch = useDispatch()
    const filePath = useSelector(selectInputNodeSelectedFilePath(nodeId))
    const onChangeFilePath = (path: string) => {
      dispatch(
        setInputImageNodeFile({
          nodeId,
          filePath: path,
        }),
      )
      const fileName = path.split('/').reverse()[0]
      dispatch(
        edifFlowElementsLabelById({
          nodeId,
          fileName,
        }),
      )
    }

    const theme = useTheme()
    const returnType = 'ImageData'
    const imageColor = useHandleColor(returnType)

    const onClickDeleteIcon = () => {
      dispatch(deleteFlowElementsById(nodeId))
    }

    return (
      <div
        style={{
          height: '100%',
          width: '250px',
          background: elementSelected
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
        <ImageFileSelect
          nodeId={nodeId}
          onChangeFilePath={onChangeFilePath}
          filePath={filePath ?? ''}
        />
        <Handle
          type="source"
          position={Position.Right}
          id={toHandleId(nodeId, 'image', returnType)}
          style={{
            ...sourceHandleStyle,
            background: imageColor,
          }}
          isValidConnection={isValidConnection}
        />
      </div>
    )
  },
)

const ImageFileSelect = React.memo<{
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
  } = useFileUploader(FILE_TYPE_SET.IMAGE)
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
        fileTreeType={FILE_TREE_TYPE_SET.IMAGE}
        selectButtonLabel="Select Image"
      />
      {error != null && (
        <Typography variant="caption" color="error">
          {error}
        </Typography>
      )}
    </>
  )
})

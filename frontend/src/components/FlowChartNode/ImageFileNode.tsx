import React, { CSSProperties, useState } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { alpha, useTheme } from '@material-ui/core'
import TextField from '@material-ui/core/TextField'
import { Handle, Position, NodeProps } from 'react-flow-renderer'

import { FlexLayoutModelContext } from 'App'
import { useTabAction } from 'FlexLayoutHook'
import { OUTPUT_TABSET_ID } from 'const/flexlayout'
import { uploadImageFile } from 'store/slice/FileData/FileDataAction'
import {
  imageIsUploadingByIdSelector,
  imageUploadingProgressSelector,
} from 'store/slice/FileData/FileDataSelector'
import { selectImageFile } from 'store/slice/FileData/FileData'
import { FileSelect } from './FileSelect'
import { FILE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'

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
  const [inputMaxIndex, setInputMaxIndex] = useState(10)
  const model = React.useContext(FlexLayoutModelContext)
  const actionForImageTab = useTabAction(element.id)
  const onClick = () => {
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
    }
  }
  const onChangeNumber = (event: React.ChangeEvent<HTMLInputElement>) => {
    setInputMaxIndex(Math.max(1, Number(event.target.value)))
  }
  const theme = useTheme()
  const imageIsUploading = useSelector(imageIsUploadingByIdSelector(element.id))
  const uploadProgress = useSelector(imageUploadingProgressSelector(element.id))
  return (
    <div
      className="imageFileNode"
      style={{
        background: element.selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
      onClick={onClick}
    >
      {imageIsUploading && uploadProgress != null && (
        <div style={{ marginLeft: 2, marginRight: 2 }}>
          <LinearProgressWithLabel value={uploadProgress} />
        </div>
      )}
      <ImageFileSelect nodeId={element.id} maxIndex={inputMaxIndex} />
      <TextField
        type="number"
        InputLabelProps={{
          shrink: true,
        }}
        value={inputMaxIndex}
        onChange={onChangeNumber}
      />
      <Handle
        type="source"
        position={Position.Right}
        id={`image-${element.id}`}
        style={sourceHandleStyle}
      />
    </div>
  )
})

const ImageFileSelect = React.memo<{ nodeId: string; maxIndex: number }>(
  ({ nodeId, maxIndex }) => {
    const dispatch = useDispatch()
    const model = React.useContext(FlexLayoutModelContext)
    const actionForImageTab = useTabAction(nodeId)
    const onUploadFile = (formData: FormData, fileName: string) => {
      dispatch(
        uploadImageFile({
          nodeId,
          fileName,
          formData,
          maxIndex,
        }),
      )
      if (actionForImageTab != null) {
        model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
      }
    }
    const onSelectFile = (path: string) => {
      dispatch(selectImageFile({ nodeId, path, maxIndex }))
      if (actionForImageTab != null) {
        model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
      }
    }
    return (
      <FileSelect
        nodeId={nodeId}
        onSelectFile={onSelectFile}
        onUploadFile={onUploadFile}
        fileType={FILE_TYPE_SET.IMAGE}
        selectButtonLabel="Select Image"
      />
    )
  },
)

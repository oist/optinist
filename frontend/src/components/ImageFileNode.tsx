import React, { CSSProperties, useState } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import {
  alpha,
  Button,
  LinearProgress,
  Typography,
  useTheme,
} from '@material-ui/core'
import ButtonGroup from '@material-ui/core/ButtonGroup'
import TextField from '@material-ui/core/TextField'
import { Handle, Position, NodeProps } from 'react-flow-renderer'

import { FlexLayoutModelContext } from 'App'
import { useTabAction } from 'FlexLayoutHook'
import { OUTPUT_TABSET_ID } from 'const/flexlayout'
import { uploadImageFile } from 'store/slice/ImageFile/ImageFileAction'
import { imageIsUploadingByIdSelector } from 'store/slice/ImageFile/ImageFileSelector'
import { filePathSelector } from 'store/slice/Element/ElementSelector'
import { selectImageFile } from 'store/slice/ImageFile/ImageFile'
import { FileSelectDialog } from 'components/FileSelectDialog'

export const ImageFileNode = React.memo<NodeProps>((element) => {
  const targetHandleStyle: CSSProperties = {
    width: 8,
    height: '100%',
    border: '1px solid',
    borderColor: '#555',
    borderRadius: 0,
  }
  const sourceHandleStyle: CSSProperties = { ...targetHandleStyle }
  const [inputFileNumber, setInputFileNumber] = useState(10)

  const model = React.useContext(FlexLayoutModelContext)
  const actionForImageTab = useTabAction(element.id)
  const onClick = () => {
    if (actionForImageTab != null) {
      model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
    }
  }

  const theme = useTheme()

  const onChangeNumber = (event: React.ChangeEvent<HTMLInputElement>) => {
    setInputFileNumber(Math.max(1, Number(event.target.value)))
  }

  const imageIsUploading = useSelector(imageIsUploadingByIdSelector(element.id))
  return (
    <div
      style={{
        height: '100%',
        background: element.selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
      onClick={onClick}
    >
      {imageIsUploading && <LinearProgress />}
      <FileSelect nodeId={element.id} maxFileNum={inputFileNumber} />
      <TextField
        type="number"
        InputLabelProps={{
          shrink: true,
        }}
        value={inputFileNumber}
        onChange={onChangeNumber}
      />
      <Handle
        type="source"
        position={Position.Right}
        id="a"
        style={sourceHandleStyle}
      />
    </div>
  )
})

const FileSelect = React.memo<{ nodeId: string; maxFileNum: number }>(
  ({ nodeId, maxFileNum }) => {
    const dispatch = useDispatch()
    const filePath = useSelector(filePathSelector(nodeId))
    const model = React.useContext(FlexLayoutModelContext)
    const actionForImageTab = useTabAction(nodeId)
    const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
      event.preventDefault()
      if (event.target.files != null && event.target.files[0] != null) {
        const file = event.target.files[0]
        const formData = new FormData()
        formData.append('file', file)
        const fileName = file.name
        dispatch(
          uploadImageFile({
            nodeId,
            fileName,
            formData,
            inputFileNumber: maxFileNum,
          }),
        )
        if (actionForImageTab != null) {
          model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
        }
      }
    }
    const inputRef = React.useRef<HTMLInputElement>(null)
    const onClick = () => {
      if (inputRef.current != null) {
        inputRef.current.click()
      }
    }
    const [open, setOpen] = React.useState(false)
    const onSelectFile = (path: string) => {
      dispatch(selectImageFile({ nodeId, path, maxIndex: maxFileNum }))
      setOpen(false)
      if (actionForImageTab != null) {
        model.doAction(actionForImageTab('image', OUTPUT_TABSET_ID))
      }
    }
    return (
      <div
        style={{
          padding: 5,
        }}
      >
        <ButtonGroup size="small" style={{ marginRight: 4 }}>
          <Button variant="outlined" onClick={() => setOpen(true)}>
            ファイルを選択
          </Button>
          <Button onClick={onClick} variant="outlined">
            またはアップロード
          </Button>
        </ButtonGroup>
        <div>
          <input
            ref={inputRef}
            type="file"
            onChange={onFileInputChange}
            accept=".tif"
            style={{
              visibility: 'hidden',
              width: 0,
              height: 0,
            }}
          />
          <Typography variant="caption">
            {filePath != null ? filePath : '---'}
          </Typography>
        </div>
        <FileSelectDialog
          selectedFilePath={filePath ?? ''}
          open={open}
          onClickOk={onSelectFile}
          onClickCancel={() => {
            setOpen(false)
          }}
          onClose={() => setOpen(false)}
          fileType="image"
        />
      </div>
    )
  },
)

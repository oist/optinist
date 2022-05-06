import React from 'react'
import { Button, Typography } from '@mui/material'
import ButtonGroup from '@mui/material/ButtonGroup'

import { FileSelectDialog } from 'components/FileSelectDialog'
import {
  FILE_TREE_TYPE,
  FILE_TREE_TYPE_SET,
} from 'store/slice/FilesTree/FilesTreeType'
import { LinearProgressWithLabel } from './LinerProgressWithLabel'
import { FILE_TYPE } from 'store/slice/InputNode/InputNodeType'
import { useFileUploader } from 'store/slice/FileUploader/FileUploaderHook'
import { getLabelByPath } from 'store/slice/FlowElement/FlowElementUtils'

export const FileSelect = React.memo<{
  multiSelect?: boolean
  filePath: string | string[]
  fileType: FILE_TYPE
  nodeId?: string
  onChangeFilePath: (path: string | string[]) => void
}>(({ multiSelect = false, filePath, nodeId, fileType, onChangeFilePath }) => {
  const {
    // filePath: uploadedFilePath,
    onUploadFile,
    pending,
    uninitialized,
    progress,
    error,
  } = useFileUploader({ fileType, nodeId })
  const onUploadFileHandle = (formData: FormData, fileName: string) => {
    onUploadFile(formData, fileName)
  }
  return (
    <>
      {!uninitialized && pending && progress != null && (
        <div style={{ marginLeft: 2, marginRight: 2 }}>
          <LinearProgressWithLabel value={progress} />
        </div>
      )}
      <FileSelectImple
        multiSelect={multiSelect}
        filePath={filePath}
        onSelectFile={onChangeFilePath}
        onUploadFile={onUploadFileHandle}
        fileTreeType={fileType}
        selectButtonLabel={`Select ${fileType}`}
      />
      {error != null && (
        <Typography variant="caption" color="error">
          {error}
        </Typography>
      )}
    </>
  )
})

type FileSelectImpleProps = {
  multiSelect?: boolean
  filePath: string | string[]
  onUploadFile: (formData: FormData, fileName: string) => void
  onSelectFile: (path: string | string[]) => void
  fileTreeType?: FILE_TREE_TYPE
  selectButtonLabel?: string
  uploadButtonLabel?: string
}

export const FileSelectImple = React.memo<FileSelectImpleProps>(
  ({
    multiSelect = false,
    filePath,
    onSelectFile,
    onUploadFile,
    fileTreeType,
    selectButtonLabel,
    uploadButtonLabel,
  }) => {
    const onFileInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
      event.preventDefault()
      if (event.target.files != null && event.target.files[0] != null) {
        const file = event.target.files[0]
        const formData = new FormData()
        formData.append('file', file)
        const fileName = file.name
        onUploadFile(formData, fileName)
      }
    }
    const inputRef = React.useRef<HTMLInputElement>(null)
    const onClick = () => {
      if (inputRef.current != null) {
        inputRef.current.click()
      }
    }
    const [open, setOpen] = React.useState(false)
    const accept = getFileInputAccept(fileTreeType)
    const fileName = getLabelByPath(filePath)
    return (
      <div
        style={{
          padding: 5,
        }}
      >
        <ButtonGroup size="small" style={{ marginRight: 4 }}>
          <Button variant="outlined" onClick={() => setOpen(true)}>
            {!!selectButtonLabel ? selectButtonLabel : 'Select File'}
          </Button>
          <Button onClick={onClick} variant="outlined">
            {!!uploadButtonLabel ? uploadButtonLabel : 'or Upload'}
          </Button>
        </ButtonGroup>
        <div>
          <input
            ref={inputRef}
            type="file"
            onChange={onFileInputChange}
            accept={accept}
            style={{
              visibility: 'hidden',
              width: 0,
              height: 0,
            }}
          />
          <Typography className="selectFilePath" variant="caption">
            {!!fileName ? fileName : 'No file is selected.'}
          </Typography>
        </div>
        <FileSelectDialog
          multiSelect={multiSelect}
          initialFilePath={filePath}
          open={open}
          onClickOk={(path) => {
            onSelectFile(path)
            setOpen(false)
          }}
          onClickCancel={() => {
            setOpen(false)
          }}
          onClose={() => setOpen(false)}
          fileType={fileTreeType}
        />
      </div>
    )
  },
)

function getFileInputAccept(fileType: FILE_TREE_TYPE | undefined) {
  switch (fileType) {
    case FILE_TREE_TYPE_SET.IMAGE:
      return '.tif'
    case FILE_TREE_TYPE_SET.CSV:
      return '.csv'
    case FILE_TREE_TYPE_SET.HDF5:
      return '.hdf5,.nwb'
    default:
      return undefined
  }
}

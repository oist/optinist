import React from 'react'
import { Button, Typography } from '@material-ui/core'
import ButtonGroup from '@material-ui/core/ButtonGroup'

import { FileSelectDialog } from 'components/FileSelectDialog'
import {
  FILE_TREE_TYPE,
  FILE_TREE_TYPE_SET,
} from 'store/slice/FilesTree/FilesTreeType'

type FileSelectProps = {
  filePath: string
  onUploadFile: (formData: FormData, fileName: string) => void
  onSelectFile: (path: string) => void
  fileTreeType?: FILE_TREE_TYPE
  selectButtonLabel?: string
  uploadButtonLabel?: string
}

export const FileSelect = React.memo<FileSelectProps>(
  ({
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
    const fileName = filePath.split('/').reverse()[0]
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
          selectedFilePath={fileName ?? ''}
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
    default:
      return undefined
  }
}

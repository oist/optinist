import React from 'react'
import { useSelector } from 'react-redux'
import { Button, Typography } from '@material-ui/core'
import ButtonGroup from '@material-ui/core/ButtonGroup'

import { filePathSelector } from 'store/slice/Element/ElementSelector'
import { FileSelectDialog } from 'components/FileSelectDialog'
import { FILE_TYPE, FILE_TYPE_SET } from 'store/slice/FilesTree/FilesTreeType'

type FileSelectProps = {
  nodeId: string
  onUploadFile: (formData: FormData, fileName: string) => void
  onSelectFile: (path: string) => void
  fileType?: FILE_TYPE
  selectButtonLabel?: string
  uploadButtonLabel?: string
}

export const FileSelect = React.memo<FileSelectProps>(
  ({
    nodeId,
    onSelectFile,
    onUploadFile,
    fileType,
    selectButtonLabel,
    uploadButtonLabel,
  }) => {
    const filePath = useSelector(filePathSelector(nodeId))
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
    const accept = getFileInputAccept(fileType)
    return (
      <div
        style={{
          padding: 5,
        }}
      >
        <ButtonGroup size="small" style={{ marginRight: 4 }}>
          <Button variant="outlined" onClick={() => setOpen(true)}>
            {!!selectButtonLabel ? selectButtonLabel : 'ファイルを選択'}
          </Button>
          <Button onClick={onClick} variant="outlined">
            {!!uploadButtonLabel ? uploadButtonLabel : 'またはアップロード'}
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
          <Typography variant="caption">
            {filePath != null ? filePath : '---'}
          </Typography>
        </div>
        <FileSelectDialog
          selectedFilePath={filePath ?? ''}
          open={open}
          onClickOk={(path) => {
            onSelectFile(path)
            setOpen(false)
          }}
          onClickCancel={() => {
            setOpen(false)
          }}
          onClose={() => setOpen(false)}
          fileType={fileType}
        />
      </div>
    )
  },
)

function getFileInputAccept(fileType: FILE_TYPE | undefined) {
  switch (fileType) {
    case FILE_TYPE_SET.IMAGE:
      return '.tif'
    case FILE_TYPE_SET.CSV:
      return '.csv'
    default:
      return undefined
  }
}

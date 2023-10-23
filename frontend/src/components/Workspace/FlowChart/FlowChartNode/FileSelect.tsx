import { ChangeEvent, memo, useContext, useRef } from "react"

import { Button, Tooltip, Typography } from "@mui/material"
import ButtonGroup from "@mui/material/ButtonGroup"

import { FILE_TREE_TYPE, FILE_TREE_TYPE_SET } from "api/files/Files"
import { DialogContext } from "components/Workspace/FlowChart/DialogContext"
import { LinearProgressWithLabel } from "components/Workspace/FlowChart/FlowChartNode/LinerProgressWithLabel"
import { useFileUploader } from "store/slice/FileUploader/FileUploaderHook"
import { getLabelByPath } from "store/slice/FlowElement/FlowElementUtils"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"

interface FileSelectProps {
  multiSelect?: boolean
  filePath: string | string[]
  fileType: FILE_TYPE
  nodeId?: string
  onChangeFilePath: (path: string | string[]) => void
}

export const FileSelect = memo(function FileSelect({
  multiSelect = false,
  filePath,
  nodeId,
  fileType,
  onChangeFilePath,
}: FileSelectProps) {
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

interface FileSelectImpleProps {
  multiSelect?: boolean
  filePath: string | string[]
  onUploadFile: (formData: FormData, fileName: string) => void
  onSelectFile: (path: string | string[]) => void
  fileTreeType?: FILE_TREE_TYPE
  selectButtonLabel?: string
  uploadButtonLabel?: string
}

export const FileSelectImple = memo(function FileSelectImple({
  multiSelect = false,
  filePath,
  onSelectFile,
  onUploadFile,
  fileTreeType,
  selectButtonLabel,
  uploadButtonLabel,
}: FileSelectImpleProps) {
  const { onOpenDialogFile } = useContext(DialogContext)

  const onFileInputChange = (event: ChangeEvent<HTMLInputElement>) => {
    event.preventDefault()
    if (event.target.files != null && event.target.files[0] != null) {
      const file = event.target.files[0]
      const formData = new FormData()
      formData.append("file", file)
      const fileName = file.name
      onUploadFile(formData, fileName)
    }
  }
  const inputRef = useRef<HTMLInputElement>(null)
  const onClick = () => {
    if (inputRef.current != null) {
      inputRef.current.click()
    }
  }
  const accept = getFileInputAccept(fileTreeType)
  const fileName = getLabelByPath(filePath)
  return (
    <div>
      <ButtonGroup size="small" style={{ marginRight: 4 }}>
        <Button
          variant="outlined"
          onClick={() => {
            onOpenDialogFile({
              open: true,
              multiSelect,
              filePath,
              fileTreeType,
              onSelectFile,
            })
          }}
        >
          {selectButtonLabel ? selectButtonLabel : "Select File"}
        </Button>
        <Button onClick={onClick} variant="outlined">
          {uploadButtonLabel ? uploadButtonLabel : "Load"}
        </Button>
      </ButtonGroup>
      <div>
        <input
          ref={inputRef}
          type="file"
          onChange={onFileInputChange}
          accept={accept}
          style={{
            visibility: "hidden",
            width: 0,
            height: 0,
          }}
        />
        <Tooltip title={fileName ? fileName : null}>
          <Typography className="selectFilePath" variant="body2">
            {fileName ? fileName : "No file is selected."}
          </Typography>
        </Tooltip>
      </div>
    </div>
  )
})

function getFileInputAccept(fileType: FILE_TREE_TYPE | undefined) {
  switch (fileType) {
    case FILE_TREE_TYPE_SET.IMAGE:
      return ".tif,.tiff"
    case FILE_TREE_TYPE_SET.CSV:
      return ".csv"
    case FILE_TREE_TYPE_SET.HDF5:
      return ".hdf5,.nwb"
    default:
      return undefined
  }
}

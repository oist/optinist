import {
  ChangeEvent,
  memo,
  ReactNode,
  useContext,
  useEffect,
  useRef,
} from "react"
import { useDispatch, useSelector } from "react-redux"

import AccountTreeIcon from "@mui/icons-material/AccountTree"
import AddLinkIcon from "@mui/icons-material/AddLink"
import AddPhotoAlternateIcon from "@mui/icons-material/AddPhotoAlternate"
import ChecklistRtlIcon from "@mui/icons-material/ChecklistRtl"
import { IconButton, Tooltip, Typography } from "@mui/material"
import ButtonGroup from "@mui/material/ButtonGroup"

import { FILE_TREE_TYPE, FILE_TREE_TYPE_SET } from "api/files/Files"
import { DialogContext } from "components/Workspace/FlowChart/Dialog/DialogContext"
import { ParamSettingDialog } from "components/Workspace/FlowChart/FlowChartNode/CsvFileNode"
import { LinearProgressWithLabel } from "components/Workspace/FlowChart/FlowChartNode/LinerProgressWithLabel"
import { getFilesTree } from "store/slice/FilesTree/FilesTreeAction"
import { useFileUploader } from "store/slice/FileUploader/FileUploaderHook"
import { getLabelByPath } from "store/slice/FlowElement/FlowElementUtils"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"
import {
  selectPipelineIsStartedSuccess,
  selectPipelineLatestUid,
} from "store/slice/Pipeline/PipelineSelectors"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

interface FileSelectProps {
  nameNode?: string
  multiSelect?: boolean
  filePath: string | string[]
  fileType: FILE_TYPE
  nodeId?: string
  onChangeFilePath: (path: string | string[]) => void
  setOpen?: (value: boolean) => void
}

export const FileSelect = memo(function FileSelect({
  nameNode,
  multiSelect = false,
  filePath,
  nodeId,
  fileType,
  setOpen,
  onChangeFilePath,
}: FileSelectProps) {
  const {
    // filePath: uploadedFilePath,
    onUploadFile,
    pending,
    uninitialized,
    progress,
    error,
    id,
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
      <Typography>{nameNode || fileType}</Typography>
      <FileSelectImple
        multiSelect={multiSelect}
        filePath={filePath}
        onSelectFile={onChangeFilePath}
        onUploadFile={onUploadFileHandle}
        fileTreeType={fileType}
        selectButtonLabel={<ChecklistRtlIcon />}
        uploadViaUrl={<AddLinkIcon />}
        nodeId={nodeId}
        id={id.current}
        setOpen={setOpen}
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
  selectButtonLabel?: ReactNode
  uploadButtonLabel?: string
  uploadViaUrl?: ReactNode
  nodeId?: string
  onUploadViaUrl?: (url: string) => void
  id: string
  setOpen?: (value: boolean) => void
}

export const FileSelectImple = memo(function FileSelectImple({
  multiSelect = false,
  filePath,
  onSelectFile,
  onUploadFile,
  fileTreeType,
  selectButtonLabel,
  uploadButtonLabel,
  uploadViaUrl,
  nodeId,
  id,
  setOpen,
}: FileSelectImpleProps) {
  const {
    onOpenFileSelectDialog,
    onOpenClearWorkflowIdDialog,
    onOpenInputUrlDialog,
  } = useContext(DialogContext)
  const dispatch = useDispatch<AppDispatch>()
  const currentWorkflowId = useSelector(selectPipelineLatestUid)
  const isPending = useSelector(selectPipelineIsStartedSuccess)
  const workspaceId = useSelector(selectCurrentWorkspaceId)

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
  const clickInput = () => {
    if (inputRef.current != null) {
      inputRef.current.click()
    }
  }

  const onClick = () => {
    if (currentWorkflowId != null) {
      onOpenClearWorkflowIdDialog({
        open: true,
        handleOk: () => {
          clickInput()
        },
        handleCancel: () => {},
      })
    } else {
      clickInput()
    }
  }

  const onClickViaUrl = () => {
    if (!nodeId) return
    onOpenInputUrlDialog({
      fileType: fileTreeType,
      open: true,
      filePath: filePath,
      nodeId,
      requestId: id,
    })
  }

  useEffect(() => {
    if (!nodeId) return
    onOpenInputUrlDialog({
      fileType: fileTreeType,
      open: false,
      filePath: filePath,
      nodeId,
      requestId: id,
    })
    //eslint-disable-next-line
  }, [])

  const accept = getFileInputAccept(fileTreeType)
  const fileName = getLabelByPath(filePath)

  return (
    <div>
      <ButtonGroup size="small" style={{ marginRight: 4 }}>
        <Tooltip title={"Select from uploaded files"}>
          <span>
            <IconButton
              color={"primary"}
              disabled={!!isPending}
              onClick={() => {
                onOpenFileSelectDialog({
                  open: true,
                  multiSelect,
                  filePath,
                  fileTreeType,
                  onSelectFile,
                })
                if (workspaceId && fileTreeType) {
                  dispatch(
                    getFilesTree({
                      workspaceId,
                      fileType: fileTreeType,
                    }),
                  )
                }
              }}
            >
              {selectButtonLabel ? selectButtonLabel : "Select File"}
            </IconButton>
          </span>
        </Tooltip>
        <Tooltip title={"Upload file"}>
          <span>
            <IconButton
              onClick={onClick}
              color={"primary"}
              disabled={!!isPending}
            >
              {uploadButtonLabel ? (
                uploadButtonLabel
              ) : (
                <AddPhotoAlternateIcon />
              )}
            </IconButton>
          </span>
        </Tooltip>
        {uploadViaUrl ? (
          <Tooltip title={"Upload file via URL"}>
            <span>
              <IconButton
                onClick={onClickViaUrl}
                color={"primary"}
                disabled={!!isPending}
              >
                {uploadViaUrl}
              </IconButton>
            </span>
          </Tooltip>
        ) : null}
        {fileTreeType === FILE_TREE_TYPE_SET.CSV && !!filePath && !!nodeId && (
          <Tooltip title={"Settings"}>
            <span>
              <IconButton disabled={!!isPending}>
                <ParamSettingDialog
                  nodeId={nodeId}
                  filePath={filePath as string}
                />
              </IconButton>
            </span>
          </Tooltip>
        )}
        {(
          [FILE_TREE_TYPE_SET.HDF5, FILE_TREE_TYPE_SET.MATLAB] as string[]
        ).includes(fileTreeType as string) &&
          !!filePath &&
          !!nodeId && (
            <Tooltip title={"Structure"}>
              <span>
                <IconButton
                  color={"primary"}
                  disabled={!!isPending}
                  onClick={() => setOpen?.(true)}
                >
                  <AccountTreeIcon />
                </IconButton>
              </span>
            </Tooltip>
          )}
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
        <Tooltip title={fileName ? fileName : null} placement="right">
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
    case FILE_TREE_TYPE_SET.MATLAB:
      return ".mat"
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

import { memo, useState } from "react"
import { DndProvider } from "react-dnd"
import { HTML5Backend } from "react-dnd-html5-backend"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar, VariantType } from "notistack"

import { Box, FormHelperText, Popover } from "@mui/material"
import { grey } from "@mui/material/colors"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { SectionTitle } from "components/common/ParamSection"
import PopupInputUrl from "components/PopupInputUrl"
import { AlgorithmOutputDialog } from "components/Workspace/FlowChart/Dialog/AlgorithmOutputDialog"
import { ClearWorkflowIdDialog } from "components/Workspace/FlowChart/Dialog/ClearWorkflowIdDialog"
import {
  ClearWorkflowIdDialogValue,
  ErrorDialogValue,
  FileSelectDialogValue,
  DialogContext,
  FileInputUrl,
} from "components/Workspace/FlowChart/Dialog/DialogContext"
import { FileSelectDialog } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import { ReactFlowComponent } from "components/Workspace/FlowChart/ReactFlowComponent"
import RightDrawer from "components/Workspace/FlowChart/RightDrawer"
import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import { CONTENT_HEIGHT, DRAWER_WIDTH, RIGHT_DRAWER_WIDTH } from "const/Layout"
import {
  getStatusLoadViaUrl,
  uploadViaUrl,
} from "store/slice/FileUploader/FileUploaderActions"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"
import { clearCurrentPipeline } from "store/slice/Pipeline/PipelineSlice"
import { selectRightDrawerIsOpen } from "store/slice/RightDrawer/RightDrawerSelectors"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

const initDialogFile = {
  filePath: "",
  open: false,
  fileTreeType: undefined,
  multiSelect: false,
  onSelectFile: () => null,
}

const initClearWorkflow = {
  open: false,
  handleOk: () => null,
  handleCancel: () => null,
}

const FlowChart = memo(function FlowChart(props: UseRunPipelineReturnType) {
  const dispatch = useDispatch<AppDispatch>()

  const open = useSelector(selectRightDrawerIsOpen)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const [dialogNodeId, setDialogNodeId] = useState("")
  const [dialogFile, setDialogFile] =
    useState<FileSelectDialogValue>(initDialogFile)
  const [dialogClearWorkflowId, setDialogClearWorkflowId] =
    useState<ClearWorkflowIdDialogValue>(initClearWorkflow)
  const [messageError, setMessageError] = useState<ErrorDialogValue>({
    anchorElRef: { current: null },
    message: "",
  })
  const [dialogViaUrl, setDialogViaUrl] = useState<FileInputUrl>({
    open: false,
    nodeId: "",
    requestId: "",
  })
  const [fileViaUrl, setFileViaUrl] = useState("")
  const [errorUrl, setErrorUrl] = useState("")

  const { enqueueSnackbar } = useSnackbar()

  const handleClickVariant = (variant: VariantType, mess: string) => {
    enqueueSnackbar(mess, { variant })
  }

  const onLoadFileViaUrl = async () => {
    if (!workspaceId || !fileViaUrl) return
    try {
      setDialogViaUrl({ ...dialogViaUrl, open: false })
      const data = await dispatch(
        uploadViaUrl({
          workspaceId,
          url: fileViaUrl,
          requestId: dialogViaUrl.requestId,
        }),
      )
      let firstCall = true
      const makeApiCall = async () => {
        const statusData = await dispatch(
          getStatusLoadViaUrl({
            workspaceId,
            file_name: (data.payload as { file_name: string }).file_name,
            requestId: dialogViaUrl.requestId,
          }),
        )
        if (
          (statusData.payload as { current: number }).current !==
          (statusData.payload as { total: number }).total
        ) {
          if (!firstCall) {
            setTimeout(makeApiCall, 1000)
          } else {
            firstCall = false
            makeApiCall()
          }
        } else {
          dispatch(
            setInputNodeFilePath({
              nodeId: dialogViaUrl.nodeId,
              filePath: [(data.payload as { file_name: string }).file_name],
            }),
          )
        }
      }
      makeApiCall()
      setDialogFile({
        ...dialogFile,
        filePath: [(data.payload as { file_name: string }).file_name],
      })
    } catch {
      handleClickVariant("error", "url does not exist")
    }
  }

  return (
    <Box display="flex">
      <DialogContext.Provider
        value={{
          onOpenOutputDialog: setDialogNodeId,
          onOpenFileSelectDialog: setDialogFile,
          onOpenClearWorkflowIdDialog: setDialogClearWorkflowId,
          onOpenInputUrlDialog: setDialogViaUrl,
          onMessageError: setMessageError,
        }}
      >
        <DndProvider backend={HTML5Backend}>
          <Box width={DRAWER_WIDTH} borderRight={1} borderColor={grey[300]}>
            <Box overflow="auto" marginRight={2}>
              <CurrentPipelineInfo />
            </Box>
            <Box overflow="auto">
              <SectionTitle>Nodes</SectionTitle>
              <AlgorithmTreeView />
            </Box>
          </Box>
          <MainContents open={open}>
            <ReactFlowComponent {...props} />
            {dialogNodeId && (
              <AlgorithmOutputDialog
                nodeId={dialogNodeId}
                open
                onClose={() => setDialogNodeId("")}
              />
            )}
            {dialogFile.open && (
              <FileSelectDialog
                multiSelect={dialogFile.multiSelect}
                initialFilePath={dialogFile.filePath}
                open={dialogFile.open}
                onClickOk={(path: string | string[]) => {
                  dialogFile.onSelectFile(path)
                  setDialogFile(initDialogFile)
                }}
                onClickCancel={() => {
                  setDialogFile(initDialogFile)
                }}
                fileType={dialogFile.fileTreeType}
              />
            )}
            {dialogClearWorkflowId.open && (
              <ClearWorkflowIdDialog
                open={dialogClearWorkflowId.open}
                onConfirm={() => {
                  dispatch(clearCurrentPipeline())
                  dialogClearWorkflowId.handleOk()
                  setDialogClearWorkflowId(initClearWorkflow)
                }}
                onCancel={() => {
                  dialogClearWorkflowId.handleCancel()
                  setDialogClearWorkflowId(initClearWorkflow)
                }}
              />
            )}
            {dialogViaUrl.fileType ? (
              <PopupInputUrl
                open={dialogViaUrl.open}
                value={fileViaUrl}
                setValue={setFileViaUrl}
                handleClose={() => {
                  setFileViaUrl("")
                  setDialogViaUrl({ ...dialogViaUrl, open: false })
                }}
                onLoadFileViaUrl={onLoadFileViaUrl}
                setError={setErrorUrl}
                error={errorUrl}
                fileType={dialogViaUrl.fileType}
              />
            ) : null}
            {messageError?.message && (
              <Popover
                open
                anchorEl={messageError.anchorElRef.current}
                onClose={() =>
                  setMessageError({
                    anchorElRef: { current: null },
                    message: "",
                  })
                }
                anchorOrigin={{
                  vertical: "top",
                  horizontal: "right",
                }}
                transformOrigin={{
                  vertical: "bottom",
                  horizontal: "left",
                }}
              >
                <div style={{ margin: 8 }}>
                  <FormHelperText error={true}>
                    {messageError.message}
                  </FormHelperText>
                </div>
              </Popover>
            )}
          </MainContents>
        </DndProvider>
        <RightDrawer />
      </DialogContext.Provider>
    </Box>
  )
})

const MainContents = styled("main")<{ open: boolean }>(
  ({ theme }) => ({
    flexDirection: "column",
    flexGrow: 1,
    minHeight: CONTENT_HEIGHT,
    transition: theme.transitions.create("margin", {
      easing: theme.transitions.easing.sharp,
      duration: theme.transitions.duration.leavingScreen,
    }),
    marginRight: -RIGHT_DRAWER_WIDTH,
  }),
  ({ open, theme }) =>
    open
      ? {
          transition: theme.transitions.create("margin", {
            easing: theme.transitions.easing.easeOut,
            duration: theme.transitions.duration.enteringScreen,
          }),
          marginRight: 0,
        }
      : undefined,
)

export default FlowChart

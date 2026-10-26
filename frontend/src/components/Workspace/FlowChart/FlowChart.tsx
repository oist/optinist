import { memo, useCallback, useState } from "react"
import { DndProvider } from "react-dnd"
import { HTML5Backend } from "react-dnd-html5-backend"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar, VariantType } from "notistack"

import CachedIcon from "@mui/icons-material/Cached"
import {
  Box,
  CircularProgress,
  FormHelperText,
  IconButton,
  Popover,
  Tooltip,
} from "@mui/material"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { DevelopmentInformation } from "components/common/DevelopmentInformation"
import { LeftSidebarContainer } from "components/common/LeftSidebarContainer"
import { SectionTitle } from "components/common/ParamSection"
import { AlgorithmOutputDialog } from "components/Workspace/FlowChart/Dialog/AlgorithmOutputDialog"
import { ClearWorkflowIdDialog } from "components/Workspace/FlowChart/Dialog/ClearWorkflowIdDialog"
import {
  ClearWorkflowIdDialogValue,
  ErrorDialogValue,
  FileSelectDialogValue,
  DialogContext,
  FileInputUrl,
  RoiSelectedProvider,
} from "components/Workspace/FlowChart/Dialog/DialogContext"
import { FileSelectDialog } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import ModalLogs from "components/Workspace/FlowChart/ModalLogs"
import { ReactFlowComponent } from "components/Workspace/FlowChart/ReactFlowComponent"
import RightDrawer from "components/Workspace/FlowChart/RightDrawer"
import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import PopupInputUrl from "components/Workspace/PopupInputUrl"
import { CONTENT_HEIGHT, RIGHT_DRAWER_WIDTH } from "const/Layout"
import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"
import {
  getStatusLoadViaUrl,
  uploadViaUrl,
} from "store/slice/FileUploader/FileUploaderActions"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import { selectLogsModalIsOpen } from "store/slice/LogsModal/LogsModalSelectors"
import {
  closeLogsModal,
  openLogsModal,
} from "store/slice/LogsModal/LogsModalSlice"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"
import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
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
  const isPending = useSelector(selectPipelineIsStartedSuccess)

  const open = useSelector(selectRightDrawerIsOpen)
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const logsModalOpen = useSelector(selectLogsModalIsOpen)
  const isDevelopment = process.env.NODE_ENV === "development"

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
  const [dialogFilterNodeId, setFilterDialogNodeId] = useState("")
  const [nodeRefresh, setNodeRefresh] = useState(false)
  const [leftSidebarOpen, setLeftSidebarOpen] = useState(true)

  const handleRefreshAlgoList = async () => {
    setNodeRefresh(true)

    try {
      await dispatch(getAlgoList()).unwrap()
      handleClickVariant("success", "Algorithm List Refreshed")
    } catch (error) {
      handleClickVariant("error", "Failed to get Algorithm List")
    } finally {
      setNodeRefresh(false)
    }
  }
  const { enqueueSnackbar } = useSnackbar()

  const handleClickVariant = (variant: VariantType, mess: string) => {
    enqueueSnackbar(mess, { variant })
  }

  const onLoadFileViaUrl = async () => {
    if (!workspaceId || !fileViaUrl) return
    setDialogViaUrl({ ...dialogViaUrl, open: false })
    const data = await dispatch(
      uploadViaUrl({
        workspaceId,
        url: fileViaUrl,
        requestId: dialogViaUrl.requestId,
      }),
    )
    if ((data as { error: { message: string } }).error) {
      handleClickVariant("error", "Url does not exist!")
      return
    }
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
        handleClickVariant("success", "File upload finished!")
        dispatch(
          setInputNodeFilePath({
            nodeId: dialogViaUrl.nodeId,
            filePath: (data.payload as { file_name: string }).file_name,
          }),
        )
      }
    }
    makeApiCall()
    setDialogFile({
      ...dialogFile,
      filePath: [(data.payload as { file_name: string }).file_name],
    })
  }

  const onCloseModalLogs = useCallback(() => {
    dispatch(closeLogsModal())
  }, [dispatch])

  const onOpenLogs = useCallback(
    (open: boolean) => {
      if (open) {
        dispatch(openLogsModal())
      } else {
        dispatch(closeLogsModal())
      }
    },
    [dispatch],
  )

  return (
    <Box display="flex" position="relative">
      <DialogContext.Provider
        value={{
          onOpenOutputDialog: setDialogNodeId,
          onOpenFileSelectDialog: setDialogFile,
          onOpenClearWorkflowIdDialog: setDialogClearWorkflowId,
          onOpenInputUrlDialog: setDialogViaUrl,
          onMessageError: setMessageError,
          onOpenFilterDialog: setFilterDialogNodeId,
          onOpenLogs,
          dialogFilterNodeId,
          isOutput: true,
        }}
      >
        <DndProvider backend={HTML5Backend}>
          <LeftSidebarContainer
            isOpen={leftSidebarOpen}
            onToggle={() => setLeftSidebarOpen(!leftSidebarOpen)}
          >
            <Box overflow="auto" marginRight={2}>
              <CurrentPipelineInfo />
            </Box>
            {isDevelopment && (
              <Box overflow="auto" marginRight={2}>
                <DevelopmentInformation />
              </Box>
            )}
            <Box overflow="auto">
              <SectionTitle
                style={{ display: "flex", alignItems: "center", gap: "8px" }}
              >
                Nodes
                <Tooltip title="Refresh Node">
                  <span>
                    <IconButton
                      onClick={handleRefreshAlgoList}
                      color="primary"
                      disabled={nodeRefresh || !!isPending} // Disable while loading
                    >
                      {nodeRefresh ? (
                        <CircularProgress size={24} />
                      ) : (
                        <CachedIcon />
                      )}
                    </IconButton>
                  </span>
                </Tooltip>
              </SectionTitle>
              <AlgorithmTreeView />
            </Box>
          </LeftSidebarContainer>
          <MainContents open={open}>
            <ReactFlowComponent {...props} />
            {dialogNodeId && (
              <AlgorithmOutputDialog
                nodeId={dialogNodeId}
                open
                onClose={() => setDialogNodeId("")}
              />
            )}
            {dialogFilterNodeId && (
              <RoiSelectedProvider>
                <AlgorithmOutputDialog
                  nodeId={dialogFilterNodeId}
                  open
                  onClose={() => setFilterDialogNodeId("")}
                />
              </RoiSelectedProvider>
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
                <div style={{ margin: 8, whiteSpace: "pre-wrap" }}>
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
      {logsModalOpen && <ModalLogs isOpen onClose={onCloseModalLogs} />}
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

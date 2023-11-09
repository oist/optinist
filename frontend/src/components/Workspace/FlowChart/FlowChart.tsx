import { memo, useState } from "react"
import { DndProvider } from "react-dnd"
import { HTML5Backend } from "react-dnd-html5-backend"
import { useDispatch, useSelector } from "react-redux"

import { Box, FormHelperText, Popover, Typography } from "@mui/material"
import { grey } from "@mui/material/colors"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { AlgorithmOutputDialog } from "components/Workspace/FlowChart/Dialog/AlgorithmOutputDialog"
import { ClearWorkflowIdDialog } from "components/Workspace/FlowChart/Dialog/ClearWorkflowIdDialog"
import {
  ClearWorkflowIdDialogValue,
  ErrorDialogValue,
  FileSelectDialogValue,
  DialogContext,
} from "components/Workspace/FlowChart/Dialog/DialogContext"
import { FileSelectDialog } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import { ReactFlowComponent } from "components/Workspace/FlowChart/ReactFlowComponent"
import RightDrawer from "components/Workspace/FlowChart/RightDrawer"
import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import { CONTENT_HEIGHT, DRAWER_WIDTH, RIGHT_DRAWER_WIDTH } from "const/Layout"
import { UseRunPipelineReturnType } from "store/slice/Pipeline/PipelineHook"
import { clearCurrentPipeline } from "store/slice/Pipeline/PipelineSlice"
import { selectRightDrawerIsOpen } from "store/slice/RightDrawer/RightDrawerSelectors"

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
  const dispatch = useDispatch()

  const open = useSelector(selectRightDrawerIsOpen)
  const [dialogNodeId, setDialogNodeId] = useState("")
  const [dialogFile, setDialogFile] =
    useState<FileSelectDialogValue>(initDialogFile)
  const [dialogClearWorkflowId, setDialogClearWorkflowId] =
    useState<ClearWorkflowIdDialogValue>(initClearWorkflow)
  const [messageError, setMessageError] = useState<ErrorDialogValue>({
    anchorElRef: { current: null },
    message: "",
  })

  return (
    <RootDiv>
      <DialogContext.Provider
        value={{
          onOpenOutputDialog: setDialogNodeId,
          onOpenFileSelectDialog: setDialogFile,
          onOpenClearWorkflowIdDialog: setDialogClearWorkflowId,
          onMessageError: setMessageError,
        }}
      >
        <DndProvider backend={HTML5Backend}>
          <Box
            sx={{
              width: DRAWER_WIDTH,
            }}
            borderRight={1}
            borderColor={grey[300]}
          >
            <CurrentPipelineInfo />
            <Typography variant="body2" color="textSecondary">
              NODES
            </Typography>
            <DrawerContents>
              <AlgorithmTreeView />
            </DrawerContents>
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
                handleOk={() => {
                  dispatch(clearCurrentPipeline())
                  dialogClearWorkflowId.handleOk()
                  setDialogClearWorkflowId(initClearWorkflow)
                }}
                handleCancel={() => {
                  dialogClearWorkflowId.handleCancel()
                  setDialogClearWorkflowId(initClearWorkflow)
                }}
              />
            )}
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
    </RootDiv>
  )
})

const RootDiv = styled("div")({
  display: "flex",
})

const DrawerContents = styled("div")({
  overflow: "auto",
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

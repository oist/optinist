import { memo, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import {
  Button,
  Dialog,
  DialogContent,
  DialogTitle,
  DialogActions,
  Switch,
  FormControlLabel,
  TextField,
  Box,
  LinearProgress,
  Typography,
} from "@mui/material"

import { FileSelect } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import { toHandleId } from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { PresentationalMatlabPlot } from "components/Workspace/Visualize/Plot/MatlabPlot"
import { HANDLE_STYLE } from "const/flowchart"
import { getMatlabData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectMatlabDataError,
  selectMatlabDataIsFulfilled,
  selectMatlabDataIsInitialized,
  selectMatlabDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  selectCsvInputNodeSelectedFilePath,
  selectInputNodeDefined,
  selectMatlabInputNodeParamSetHeader,
  selectMatlabInputNodeParamSetIndex,
  selectMatlabInputNodeParamTranspose,
} from "store/slice/InputNode/InputNodeSelectors"
import { setMatlabInputNodeParam } from "store/slice/InputNode/InputNodeSlice"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const MatlabFileNode = memo(function MatlabFileNode(element: NodeProps) {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <MatlabFileNodeImple {...element} />
  } else {
    return null
  }
})

const MatlabFileNodeImple = memo(function MatlabFileNodeImple({
  id: nodeId,
  selected,
}: NodeProps) {
  const dispatch = useDispatch<AppDispatch>()
  const filePath = useSelector(selectCsvInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowNodeById(nodeId))
  }

  return (
    <NodeContainer nodeId={nodeId} selected={selected}>
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: "black", position: "absolute", top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        onChangeFilePath={(path: string | string[]) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        fileType={FILE_TYPE_SET.CSV}
        filePath={filePath ?? ""}
      />
      {!!filePath && <ParamSettingDialog nodeId={nodeId} filePath={filePath} />}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, "matlab", "MatlabData")}
        style={{ ...HANDLE_STYLE }}
      />
    </NodeContainer>
  )
})

interface ParamSettingDialogProps extends NodeIdProps {
  filePath: string
}

export const ParamSettingDialog = memo(function ParamSettingDialog({
  nodeId,
  filePath,
}: ParamSettingDialogProps) {
  const [open, setOpen] = useState(false)
  // OK時のみStoreに反映させるため一時的な値をuseStateで保持しておく。
  // useStateの初期値はselectorで取得。
  const [setHeader, setSetHeader] = useState(
    useSelector(selectMatlabInputNodeParamSetHeader(nodeId)),
  )
  const [setIndex, setSetIndex] = useState(
    useSelector(selectMatlabInputNodeParamSetIndex(nodeId)),
  )
  const [transpose, setTranspose] = useState(
    useSelector(selectMatlabInputNodeParamTranspose(nodeId)),
  )
  const dispatch = useDispatch<AppDispatch>()
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(
      setMatlabInputNodeParam({
        nodeId,
        param: { setHeader, setIndex, transpose },
      }),
    )
  }

  return (
    <>
      <Button onClick={() => setOpen(true)} sx={{ padding: 0 }}>
        Settings
      </Button>
      <Dialog open={open}>
        <DialogTitle>Matlab Setting</DialogTitle>
        <DialogContent dividers>
          <Box sx={{ display: "flex", p: 1, m: 1, alignItems: "flex-start" }}>
            <FormControlLabel
              sx={{ margin: (theme) => theme.spacing(0, 1, 0, 1) }}
              control={
                <Switch
                  checked={transpose}
                  onChange={(event) => setTranspose(event.target.checked)}
                />
              }
              label="Transpose"
            />
            <TextField
              label="header"
              sx={{
                width: 100,
                margin: (theme) => theme.spacing(0, 1, 0, 1),
              }}
              type="number"
              InputLabelProps={{
                shrink: true,
              }}
              onChange={(event) => {
                const value = Number(event.target.value)
                if (value >= 0) {
                  setSetHeader(value)
                }
              }}
              value={setHeader}
            />
            <FormControlLabel
              sx={{ margin: (theme) => theme.spacing(0, 1, 0, 1) }}
              control={
                <Switch
                  checked={setIndex}
                  onChange={(event) => setSetIndex(event.target.checked)}
                />
              }
              label="Set Index"
            />
          </Box>
          <Typography variant="h6">Preview</Typography>
          <MatlabPreview
            filePath={filePath}
            transpose={transpose}
            setHeader={setHeader}
            setIndex={setIndex}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            cancel
          </Button>
          <Button onClick={onClickOk} color="primary" variant="outlined">
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

interface MatlabPreviewProps {
  filePath: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}

const MatlabPreview = memo(function MatlabPreview({
  filePath: path,
  ...otherProps
}: MatlabPreviewProps) {
  const isInitialized = useSelector(selectMatlabDataIsInitialized(path))
  const isPending = useSelector(selectMatlabDataIsPending(path))
  const isFulfilled = useSelector(selectMatlabDataIsFulfilled(path))
  const error = useSelector(selectMatlabDataError(path))
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  useEffect(() => {
    if (workspaceId && !isInitialized) {
      dispatch(getMatlabData({ path, workspaceId }))
    }
  }, [dispatch, isInitialized, path, workspaceId])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <PresentationalMatlabPlot path={path} {...otherProps} />
  } else {
    return null
  }
})

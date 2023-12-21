import { memo, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { Handle, Position, NodeProps } from "reactflow"

import SettingsIcon from "@mui/icons-material/Settings"
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
  IconButton,
} from "@mui/material"

import { FileSelect } from "components/Workspace/FlowChart/FlowChartNode/FileSelect"
import { toHandleId } from "components/Workspace/FlowChart/FlowChartNode/FlowChartUtils"
import { NodeContainer } from "components/Workspace/FlowChart/FlowChartNode/NodeContainer"
import { PresentationalCsvPlot } from "components/Workspace/Visualize/Plot/CsvPlot"
import { HANDLE_STYLE } from "const/flowchart"
import { getCsvData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectCsvDataError,
  selectCsvDataIsFulfilled,
  selectCsvDataIsInitialized,
  selectCsvDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import { deleteFlowNodeById } from "store/slice/FlowElement/FlowElementSlice"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import { setInputNodeFilePath } from "store/slice/InputNode/InputNodeActions"
import {
  selectCsvInputNodeParamSetHeader,
  selectCsvInputNodeParamSetIndex,
  selectCsvInputNodeParamTranspose,
  selectCsvInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from "store/slice/InputNode/InputNodeSelectors"
import { setCsvInputNodeParam } from "store/slice/InputNode/InputNodeSlice"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const CsvFileNode = memo(function CsvFileNode(element: NodeProps) {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <CsvFileNodeImple {...element} />
  } else {
    return null
  }
})

const CsvFileNodeImple = memo(function CsvFileNodeImple({
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
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, "csv", "CsvData")}
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
    useSelector(selectCsvInputNodeParamSetHeader(nodeId)),
  )
  const [setIndex, setSetIndex] = useState(
    useSelector(selectCsvInputNodeParamSetIndex(nodeId)),
  )
  const [transpose, setTranspose] = useState(
    useSelector(selectCsvInputNodeParamTranspose(nodeId)),
  )
  const dispatch = useDispatch<AppDispatch>()
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(
      setCsvInputNodeParam({
        nodeId,
        param: { setHeader, setIndex, transpose },
      }),
    )
  }

  return (
    <>
      <IconButton
        onClick={() => setOpen(true)}
        sx={{ padding: 0 }}
        color={"primary"}
      >
        <SettingsIcon />
      </IconButton>
      <Dialog open={open} onClose={onClickCancel}>
        <DialogTitle>Csv Setting</DialogTitle>
        <DialogContent dividers>
          <Box
            sx={{
              display: "flex",
              alignItems: "flex-start",
            }}
          >
            <FormControlLabel
              sx={{
                margin: (theme) => theme.spacing(0, 1, 0, 1),
                whiteSpace: "nowrap",
              }}
              control={
                <Switch
                  checked={transpose}
                  onChange={(event) => setTranspose(event.target.checked)}
                />
              }
              label="Transpose"
            />
            <FormControlLabel
              sx={{
                margin: (theme) => theme.spacing(0, 1, 0, 1),
                whiteSpace: "nowrap",
              }}
              control={
                <Switch
                  checked={setHeader != null}
                  onChange={(event) => {
                    if (event.target.checked) {
                      setSetHeader(0)
                    } else {
                      setSetHeader(null)
                    }
                  }}
                />
              }
              label="Set Header"
            />
            {setHeader != null && (
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
            )}
            <FormControlLabel
              sx={{
                margin: (theme) => theme.spacing(0, 1, 0, 1),
                whiteSpace: "nowrap",
              }}
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
          <CsvPreview
            filePath={filePath}
            transpose={transpose}
            setHeader={setHeader}
            setIndex={setIndex}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined">
            cancel
          </Button>
          <Button onClick={onClickOk} variant="contained">
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

interface CsvPreviewProps {
  filePath: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}

const CsvPreview = memo(function CsvPreview({
  filePath: path,
  ...otherProps
}: CsvPreviewProps) {
  const isInitialized = useSelector(selectCsvDataIsInitialized(path))
  const isPending = useSelector(selectCsvDataIsPending(path))
  const isFulfilled = useSelector(selectCsvDataIsFulfilled(path))
  const error = useSelector(selectCsvDataError(path))
  const dispatch = useDispatch<AppDispatch>()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  useEffect(() => {
    if (workspaceId && !isInitialized) {
      dispatch(getCsvData({ path, workspaceId }))
    }
  }, [dispatch, isInitialized, path, workspaceId])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <PresentationalCsvPlot path={path} {...otherProps} />
  } else {
    return null
  }
})

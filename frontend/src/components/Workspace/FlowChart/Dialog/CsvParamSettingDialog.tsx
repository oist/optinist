import { memo, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

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

import { PresentationalCsvPlot } from "components/Workspace/Visualize/Plot/CsvPlot"
import { getCsvData } from "store/slice/DisplayData/DisplayDataActions"
import {
  selectCsvDataError,
  selectCsvDataIsFulfilled,
  selectCsvDataIsInitialized,
  selectCsvDataIsPending,
} from "store/slice/DisplayData/DisplayDataSelectors"
import { NodeIdProps } from "store/slice/FlowElement/FlowElementType"
import {
  selectCsvInputNodeParamSetHeader,
  selectCsvInputNodeParamSetIndex,
  selectCsvInputNodeParamTranspose,
} from "store/slice/InputNode/InputNodeSelectors"
import { setCsvInputNodeParam } from "store/slice/InputNode/InputNodeSlice"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

interface CsvParamSettingDialogProps extends NodeIdProps {
  filePath: string
  disabled?: boolean
}

export const CsvParamSettingDialog = memo(function CsvParamSettingDialog({
  nodeId,
  filePath,
  disabled = false,
}: CsvParamSettingDialogProps) {
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
        color={"primary"}
        disabled={disabled}
        size="small"
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
                label="header index"
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

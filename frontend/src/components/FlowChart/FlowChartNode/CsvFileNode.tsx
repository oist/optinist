import React, { CSSProperties } from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { Handle, Position, NodeProps } from 'react-flow-renderer'
import { alpha, useTheme } from '@mui/material/styles'
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
} from '@mui/material'

import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  selectCsvInputNodeParamSetHeader,
  selectCsvInputNodeParamSetIndex,
  selectCsvInputNodeParamTranspose,
  selectCsvInputNodeSelectedFilePath,
  selectInputNodeDefined,
} from 'store/slice/InputNode/InputNodeSelectors'
import { setCsvInputNodeParam } from 'store/slice/InputNode/InputNodeSlice'
import { setInputNodeFilePath } from 'store/slice/InputNode/InputNodeActions'
import { toHandleId } from './FlowChartUtils'
import { FileSelect } from './FileSelect'
import { deleteFlowElementsById } from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectCsvDataError,
  selectCsvDataIsFulfilled,
  selectCsvDataIsInitialized,
  selectCsvDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getCsvData } from 'store/slice/DisplayData/DisplayDataActions'
import { PresentationalCsvPlot } from 'components/Visualize/Plot/CsvPlot'

const sourceHandleStyle: CSSProperties = {
  width: 8,
  height: 15,
  top: 15,
  border: '1px solid',
  borderColor: '#555',
  borderRadius: 0,
}

export const CsvFileNode = React.memo<NodeProps>((element) => {
  const defined = useSelector(selectInputNodeDefined(element.id))
  if (defined) {
    return <CsvFileNodeImple {...element} />
  } else {
    return null
  }
})

const CsvFileNodeImple = React.memo<NodeProps>(({ id: nodeId, selected }) => {
  const dispatch = useDispatch()
  const filePath = useSelector(selectCsvInputNodeSelectedFilePath(nodeId))
  const onChangeFilePath = (path: string) => {
    dispatch(setInputNodeFilePath({ nodeId, filePath: path }))
  }
  const theme = useTheme()

  const onClickDeleteIcon = () => {
    dispatch(deleteFlowElementsById(nodeId))
  }

  return (
    <div
      style={{
        height: '100%',
        width: '230px',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
    >
      <button
        className="flowbutton"
        onClick={onClickDeleteIcon}
        style={{ color: 'black', position: 'absolute', top: -10, right: 10 }}
      >
        ×
      </button>
      <FileSelect
        nodeId={nodeId}
        onChangeFilePath={(path) => {
          if (!Array.isArray(path)) {
            onChangeFilePath(path)
          }
        }}
        fileType={FILE_TYPE_SET.CSV}
        filePath={filePath ?? ''}
      />
      {!!filePath && <ParamSettingDialog nodeId={nodeId} filePath={filePath} />}
      <Handle
        type="source"
        position={Position.Right}
        id={toHandleId(nodeId, 'csv', 'CsvData')}
        style={sourceHandleStyle}
      />
    </div>
  )
})

export const ParamSettingDialog = React.memo<{
  nodeId: string
  filePath: string
}>(({ nodeId, filePath }) => {
  const [open, setOpen] = React.useState(false)
  // OK時のみStoreに反映させるため一時的な値をuseStateで保持しておく。
  // useStateの初期値はselectorで取得。
  const [setHeader, setSetHeader] = React.useState(
    useSelector(selectCsvInputNodeParamSetHeader(nodeId)),
  )
  const [setIndex, setSetIndex] = React.useState(
    useSelector(selectCsvInputNodeParamSetIndex(nodeId)),
  )
  const [transpose, setTranspose] = React.useState(
    useSelector(selectCsvInputNodeParamTranspose(nodeId)),
  )
  const dispatch = useDispatch()
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
      <Button onClick={() => setOpen(true)}>Settings</Button>
      <Dialog open={open}>
        <DialogTitle>Csv Setting</DialogTitle>
        <DialogContent dividers>
          <Box sx={{ display: 'flex', p: 1, m: 1, alignItems: 'flex-start' }}>
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
          <CsvPreview
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

const CsvPreview = React.memo<{
  filePath: string
  transpose: boolean
  setHeader: number | null
  setIndex: boolean
}>(({ filePath: path, ...otherProps }) => {
  const isInitialized = useSelector(selectCsvDataIsInitialized(path))
  const isPending = useSelector(selectCsvDataIsPending(path))
  const isFulfilled = useSelector(selectCsvDataIsFulfilled(path))
  const error = useSelector(selectCsvDataError(path))
  const dispatch = useDispatch()
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getCsvData({ path }))
    }
  }, [dispatch, isInitialized, path])
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

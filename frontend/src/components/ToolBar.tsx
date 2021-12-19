import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import {
  nodeDataListForRunSelector,
  pathIsUndefinedSelector,
} from 'store/slice/Element/ElementSelector'
import { Box, IconButton, LinearProgress } from '@material-ui/core'
import Close from '@material-ui/icons/Close'
import { reflectRunPipelineResult } from 'store/slice/Algorithm/AlgorithmAction'
import { useLazyRunPipelineQuery, RunPipeLineNodeDataType } from 'api/Run/Run'
import { AlgoNodeData, CsvNodeData, ImageNodeData } from 'const/NodeData'
import { nanoid } from '@reduxjs/toolkit'
import { SnackbarProvider, SnackbarKey, useSnackbar } from 'notistack'
import { NWB } from './NWB'

export const ToolBar = React.memo(() => (
  <SnackbarProvider
    maxSnack={5}
    action={(snackbarKey) => <SnackbarCloseButton snackbarKey={snackbarKey} />}
  >
    <ToolBarImple />
  </SnackbarProvider>
))

const SnackbarCloseButton: React.FC<{ snackbarKey: SnackbarKey }> = ({
  snackbarKey,
}) => {
  const { closeSnackbar } = useSnackbar()
  return (
    <IconButton onClick={() => closeSnackbar(snackbarKey)}>
      <Close style={{ color: 'white' }} />
    </IconButton>
  )
}

export const ToolBarImple = React.memo(() => {
  const dispatch = useDispatch()
  const { enqueueSnackbar, closeSnackbar } = useSnackbar()
  const pathIsUndefined = useSelector(pathIsUndefinedSelector)
  const nodeDataListForRun = useSelector(
    nodeDataListForRunSelector,
    nodeDataListForRunEqualityFn,
  )
  const [triggerRunPipeline, result] = useLazyRunPipelineQuery()
  const [isReady, setIsReady] = React.useState(false)
  const onRunBtnClick = () => {
    if (pathIsUndefined) {
      enqueueSnackbar('failed to read file path.', { variant: 'error' })
    } else {
      triggerRunPipeline({ nodeDataListForRun, requestId: nanoid() })
      closeSnackbar()
      setIsReady(true)
    }
  }
  React.useEffect(() => {
    if (result.data != null && !result.isFetching) {
      dispatch(
        reflectRunPipelineResult({
          dto: result.data.outputPaths ?? {},
          error:
            result.data.name != null
              ? { name: result.data.name, message: result.data.message }
              : undefined,
        }),
      )
      if (result.data.status === 'ready') {
        setIsReady(true)
      } else {
        setIsReady(false)
        enqueueSnackbar(result.data.message, {
          variant:
            result.data.status === 'error'
              ? 'error'
              : result.data.status === 'success' ||
                result.data.status === 'completed'
              ? 'success'
              : undefined,
        })
      }
    }
  }, [result, enqueueSnackbar, closeSnackbar, dispatch])
  return (
    <div style={{ width: '100%' }}>
      <Box
        display="flex"
        justifyContent="flex-end"
        style={{ paddingBottom: 4 }}
      >
        <NWB />
        <Box>
          <Button
            className="ctrl_btn"
            variant="contained"
            color="primary"
            endIcon={<PlayArrowIcon />}
            onClick={onRunBtnClick}
            disabled={isReady}
          >
            run
          </Button>
        </Box>
      </Box>
      {isReady ? <LinearProgress /> : <div style={{ height: 4 }} />}
    </div>
  )
})

function nodeDataListForRunEqualityFn(
  a: (RunPipeLineNodeDataType | undefined)[],
  b: (RunPipeLineNodeDataType | undefined)[],
) {
  return (
    a === b ||
    (a.length === b.length &&
      a.every((v, i) => nodeDataForRunEqualityFn(v, b[i])))
  )
}

function nodeDataForRunEqualityFn(
  a: RunPipeLineNodeDataType | undefined,
  b: RunPipeLineNodeDataType | undefined,
) {
  if (a !== undefined && b !== undefined) {
    return (
      (a.type === 'algo' &&
        b.type === 'algo' &&
        algoNodeDataEqualityFn(a, b)) ||
      (a.type === 'image' &&
        b.type === 'image' &&
        inputNodeDataEqualityFn(a, b))
    )
  } else {
    return a === undefined && b === undefined
  }
}

function algoNodeDataEqualityFn(a: AlgoNodeData, b: AlgoNodeData) {
  return a.label === b.label && algoNodeDataParamaEqualityFn(a.param, b.param)
}

function algoNodeDataParamaEqualityFn(
  a: AlgoNodeData['param'],
  b: AlgoNodeData['param'],
) {
  if (a !== undefined && b !== undefined) {
    const aArray = Object.entries(a)
    const bArray = Object.entries(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every(([aKey, aValue], i) => {
          const [bKey, bValue] = bArray[i]
          return bKey === aKey && bValue === aValue
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}

function inputNodeDataEqualityFn(
  a: ImageNodeData | CsvNodeData,
  b: ImageNodeData | CsvNodeData,
) {
  return a.path === b.path && a.label === b.label
}

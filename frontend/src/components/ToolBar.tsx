import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import {
  elementListForRunSelector,
  pathIsUndefinedSelector,
} from 'store/slice/Element/ElementSelector'
import { Box, IconButton, LinearProgress } from '@material-ui/core'
import Close from '@material-ui/icons/Close'
import { reflectRunPipelineResult } from 'store/slice/Algorithm/AlgorithmAction'
import { useLazyRunPipelineQuery } from 'api/Run/Run'
import { nanoid } from '@reduxjs/toolkit'
import { SnackbarProvider, SnackbarKey, useSnackbar } from 'notistack'

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
  const elementListForRun = useSelector(elementListForRunSelector)
  const [triggerRunPipeline, result] = useLazyRunPipelineQuery()
  const [isReady, setIsReady] = React.useState(false)
  const onRunBtnClick = () => {
    if (pathIsUndefined) {
      enqueueSnackbar('failed to read file path.', { variant: 'error' })
    } else if (elementListForRun.edgeList.length === 0) {
      enqueueSnackbar('there are no edges.', { variant: 'error' })
    } else {
      triggerRunPipeline({ elementListForRun, requestId: nanoid() })
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

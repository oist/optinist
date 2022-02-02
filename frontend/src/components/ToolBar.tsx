import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { nanoid } from '@reduxjs/toolkit'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import { Box, IconButton, LinearProgress } from '@material-ui/core'
import Close from '@material-ui/icons/Close'
import { SnackbarProvider, SnackbarKey, useSnackbar } from 'notistack'

import { NWBSettingButton } from './FlowChart/NWB'
import { selectNwbList } from 'store/slice/NWB/NWBSelectors'
import { selectFilePathIsUndefined } from 'store/slice/InputNode/InputNodeSelectors'
import { selectElementListForRun } from 'store/slice/FlowElement/FlowElementSelectors'
import { reflectResult } from 'store/slice/RunPipelineResult/RunPipelineResultSlice'
import { RunPipeLineContext } from './RunContext'
import { SnakemakeButton } from './FlowChart/Snakemake'

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

const ToolBarImple = React.memo(() => {
  const dispatch = useDispatch()
  const { enqueueSnackbar, closeSnackbar } = useSnackbar()
  const pathIsUndefined = useSelector(selectFilePathIsUndefined)
  const nwbParam = useSelector(selectNwbList)
  const elementListForRun = useSelector(selectElementListForRun)
  // const [triggerRunPipeline, result] = useLazyRunPipelineQuery()
  const { runPipeLine, result } = React.useContext(RunPipeLineContext)
  const [isReady, setIsReady] = React.useState(false)
  const onRunBtnClick = () => {
    if (pathIsUndefined) {
      enqueueSnackbar('failed to read file path.', { variant: 'error' })
    } else if (elementListForRun.edgeList.length === 0) {
      enqueueSnackbar('there are no edges.', { variant: 'error' })
    } else {
      // triggerRunPipeline({ elementListForRun, requestId: nanoid(), nwbParam })
      runPipeLine({ elementListForRun, requestId: nanoid(), nwbParam })
      closeSnackbar()
      setIsReady(true)
    }
  }
  React.useEffect(() => {
    if (result.data != null && !result.isFetching) {
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
        dispatch(reflectResult(result.data))
      }
    }
  }, [result, enqueueSnackbar, closeSnackbar, dispatch])
  return (
    <div style={{ width: '100%' }}>
      <Box display="flex" justifyContent="flex-end" style={{ padding: 4 }}>
        <SnakemakeButton />
        <NWBSettingButton />
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

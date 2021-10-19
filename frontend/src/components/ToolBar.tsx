import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import StopIcon from '@material-ui/icons/Stop'
import Snackbar from '@material-ui/core/Snackbar'
import CloseIcon from '@material-ui/icons/Close'

import {
  runMassageSelector,
  runStatusSelector,
} from 'redux/slice/Element/ElementSelector'
import { Box, IconButton, LinearProgress } from '@material-ui/core'
import { runPipeline, stopPipeline } from 'redux/slice/Element/ElementAction'
import { RootState } from 'redux/store'
import { RUN_STATUS } from 'redux/slice/Element/ElementType'

export const ToolBar = React.memo(() => {
  const dispatch = useDispatch()
  const isRunning = useSelector(
    (state: RootState) => runStatusSelector(state) === RUN_STATUS.RUNNING,
  )
  const onRunBtnClick = () => {
    dispatch(runPipeline())
  }
  const onStopBtnClick = () => {
    dispatch(stopPipeline())
  }
  const [dialog, setDialog] = React.useState<{
    open: boolean
    message: string
  }>({ open: false, message: '' })
  const handleClose = () => {
    setDialog({ open: false, message: '' })
  }
  const runStatus = useSelector(runStatusSelector)
  const runMessage = useSelector(runMassageSelector)
  React.useEffect(() => {
    if (runStatus === RUN_STATUS.SUCCESS || runStatus === RUN_STATUS.FAILED) {
      setDialog({ open: true, message: runMessage ?? '' })
    }
  }, [runStatus, runMessage])
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
            disabled={isRunning}
          >
            run
          </Button>
        </Box>
        <Box>
          <Button
            className="ctrl_btn"
            variant="contained"
            color="secondary"
            endIcon={<StopIcon />}
            onClick={onStopBtnClick}
          >
            stop
          </Button>
        </Box>
      </Box>
      {isRunning ? <LinearProgress /> : <div style={{ height: 4 }} />}
      <Snackbar
        autoHideDuration={5000}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'left',
        }}
        open={dialog.open}
        onClose={handleClose}
        message={dialog.message}
        action={
          <IconButton onClick={handleClose} color="inherit" size="small">
            <CloseIcon />
          </IconButton>
        }
      />
    </div>
  )
})

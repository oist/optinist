import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import Snackbar from '@material-ui/core/Snackbar'
import CloseIcon from '@material-ui/icons/Close'
import Alert from '@material-ui/lab/Alert'
import {
  runMassageSelector,
  runStatusSelector,
} from 'redux/slice/Element/ElementSelector'
import { Box, IconButton, LinearProgress } from '@material-ui/core'
import { runPipeline } from 'redux/slice/Element/ElementAction'
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
  const [snackbar, setDialog] = React.useState<{
    open: boolean
    message?: string
    severity?: 'success' | 'error'
  }>({ open: false })
  const handleClose = () => {
    setDialog((prev) => ({ ...prev, open: false }))
  }
  const runStatus = useSelector(runStatusSelector)
  const runMessage = useSelector(runMassageSelector)
  React.useEffect(() => {
    if (runStatus === RUN_STATUS.SUCCESS) {
      setDialog({ open: true, message: runMessage, severity: 'success' })
    } else if (runStatus === RUN_STATUS.FAILED) {
      setDialog({ open: true, message: runMessage, severity: 'error' })
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
      </Box>
      {isRunning ? <LinearProgress /> : <div style={{ height: 4 }} />}
      <Snackbar
        autoHideDuration={5000}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'left',
        }}
        open={snackbar.open}
        onClose={handleClose}
        action={
          <IconButton onClick={handleClose} color="inherit" size="small">
            <CloseIcon />
          </IconButton>
        }
      >
        <Alert onClose={handleClose} severity={snackbar.severity}>
          {snackbar.message}
        </Alert>
      </Snackbar>
    </div>
  )
})

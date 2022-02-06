import React from 'react'

import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import { createStyles, makeStyles, Theme } from '@material-ui/core/styles'
import CloseIcon from '@material-ui/icons/Close'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { RUN_STATUS } from 'store/slice/Pipeline/PipelineType'
import { useSnackbar } from 'notistack'

export const RunButtons = React.memo<UseRunPipelineReturnType>((props) => {
  const { status, handleCancelPipeline, handleRunPipeline } = props
  const { enqueueSnackbar } = useSnackbar()
  const onClickRun = () => {
    handleRunPipeline()
    if (status === RUN_STATUS.START_SUCCESS) {
      enqueueSnackbar('started')
    }
  }
  const onClickCancel = () => {
    handleCancelPipeline()
  }

  React.useEffect(() => {
    if (status === RUN_STATUS.FINISHED) {
      enqueueSnackbar('Finished', { variant: 'success' })
    } else if (status === RUN_STATUS.ABORTED) {
      enqueueSnackbar('Aborted', { variant: 'error' })
    } else if (status === RUN_STATUS.CANCELED) {
      enqueueSnackbar('Canceled', { variant: 'info' })
    }
  }, [status, enqueueSnackbar])

  const classes = useStyles()
  return (
    <>
      <Button
        variant="contained"
        color="primary"
        endIcon={<PlayArrowIcon />}
        className={classes.button}
        onClick={onClickRun}
        disabled={status === RUN_STATUS.START_SUCCESS}
      >
        Run
      </Button>
      <Button
        variant="outlined"
        endIcon={<CloseIcon />}
        className={classes.button}
        onClick={onClickCancel}
      >
        Cancel
      </Button>
    </>
  )
})

const useStyles = makeStyles((theme: Theme) =>
  createStyles({
    button: {
      margin: theme.spacing(1),
    },
  }),
)

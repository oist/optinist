import React from 'react'

import Button from '@mui/material/Button'
import PlayArrowIcon from '@mui/icons-material/PlayArrow'
import CloseIcon from '@mui/icons-material/Close'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { RUN_STATUS } from 'store/slice/Pipeline/PipelineType'
import { useSnackbar } from 'notistack'

export const RunButtons = React.memo<UseRunPipelineReturnType>((props) => {
  const {
    status,
    filePathIsUndefined,
    handleCancelPipeline,
    handleRunPipeline,
  } = props
  const { enqueueSnackbar } = useSnackbar()
  const onClickRun = () => {
    if (!filePathIsUndefined) {
      handleRunPipeline()
    } else {
      enqueueSnackbar('please select input file', { variant: 'error' })
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

  return (
    <>
      <Button
        variant="contained"
        color="primary"
        endIcon={<PlayArrowIcon />}
        onClick={onClickRun}
        disabled={status === RUN_STATUS.START_SUCCESS}
        sx={{
          margin: (theme) => theme.spacing(1),
        }}
      >
        Run
      </Button>
      <Button
        variant="outlined"
        endIcon={<CloseIcon />}
        onClick={onClickCancel}
        sx={{
          margin: (theme) => theme.spacing(1),
        }}
      >
        Cancel
      </Button>
    </>
  )
})

import React from 'react'

import { useDispatch, useSelector } from 'react-redux'
import Button from '@material-ui/core/Button'
import PlayArrowIcon from '@material-ui/icons/PlayArrow'
import { createStyles, makeStyles, Theme } from '@material-ui/core/styles'
import CloseIcon from '@material-ui/icons/Close'

import { run } from 'store/slice/Pipeline/PilepineActions'

export const RunButtons: React.FC = () => {
  const dispatch = useDispatch()
  const onClickRun = () => {
    dispatch(run())
  }
  const onClickCancel = () => {
    //  todo
  }
  const classes = useStyles()
  return (
    <>
      <Button
        variant="contained"
        color="secondary"
        endIcon={<PlayArrowIcon />}
        className={classes.button}
        onClick={onClickRun}
      >
        Run(new)
      </Button>
      <Button
        variant="contained"
        endIcon={<CloseIcon />}
        className={classes.button}
        onClick={onClickCancel}
      >
        Cancel
      </Button>
    </>
  )
}

const useStyles = makeStyles((theme: Theme) =>
  createStyles({
    button: {
      margin: theme.spacing(1),
    },
  }),
)

import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@mui/material/Button'
import IconButton from '@mui/material/IconButton'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogTitle from '@mui/material/DialogTitle'
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline'
import { selectExperimentName } from 'store/slice/Experiments/ExperimentsSelectors'
import { deleteExperimentByUid } from 'store/slice/Experiments/ExperimentsActions'
import { ExperimentUidContext } from '../ExperimentTable'

export const DeleteButton = React.memo(() => {
  const dispatch = useDispatch()
  const uid = React.useContext(ExperimentUidContext)

  const name = useSelector(selectExperimentName(uid))
  const [open, setOpen] = React.useState(false)

  const onClickOpen = () => {
    setOpen(true)
  }
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(deleteExperimentByUid(uid))
  }
  return (
    <>
      <IconButton onClick={onClickOpen}>
        <DeleteOutlineIcon color="error" />
      </IconButton>
      <Dialog open={open}>
        <DialogTitle>Are you sure you want to delete {name}?</DialogTitle>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            Cancel
          </Button>
          <Button onClick={onClickOk} variant="outlined" autoFocus>
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  Tooltip,
} from '@mui/material'
import ClearIcon from '@mui/icons-material/Clear'
import React, { SetStateAction } from 'react'
import { useDispatch } from 'react-redux'
import { clearFlowElements } from 'store/slice/FlowElement/FlowElementSlice'
import DeleteIcon from '@mui/icons-material/Delete'

export const ClearWorkflowButton = React.memo(() => {
  const [open, setOpen] = React.useState(false)
  const openDialog = () => {
    setOpen(true)
  }

  return (
    <>
      <Tooltip title="Clear workflow">
        <IconButton onClick={openDialog}>
          <DeleteIcon color="primary" />
        </IconButton>
      </Tooltip>
      <ConfirmClearDialog open={open} setOpen={setOpen} />
    </>
  )
})

const ConfirmClearDialog = React.memo<{
  open: boolean
  setOpen: React.Dispatch<SetStateAction<boolean>>
}>(({ open, setOpen }) => {
  const dispatch = useDispatch()
  const handleClose = () => {
    setOpen(false)
  }
  const handleClear = () => {
    dispatch(clearFlowElements())
    setOpen(false)
  }

  return (
    <Dialog open={open} onClose={handleClose}>
      <DialogTitle>Confirm clear workflow</DialogTitle>
      <DialogContent>
        Are you sure you want to clear the current workflow?
      </DialogContent>
      <DialogActions>
        <Button variant="outlined" onClick={handleClose}>Cancel</Button>
        <Button variant="contained" onClick={handleClear}>Clear</Button>
      </DialogActions>
    </Dialog>
  )
})

import { Dispatch, memo, SetStateAction, useState } from "react"
import { useDispatch } from "react-redux"

import { AddToPhotos } from "@mui/icons-material"
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  IconButton,
  Tooltip,
} from "@mui/material"

import { clearFlowElements } from "store/slice/FlowElement/FlowElementSlice"

export const CreateWorkflowButton = memo(function CreateWorkflowButton() {
  const [open, setOpen] = useState(false)
  const openDialog = () => {
    setOpen(true)
  }

  return (
    <>
      <Tooltip title="Create new workflow">
        <IconButton onClick={openDialog}>
          <AddToPhotos color="primary" />
        </IconButton>
      </Tooltip>
      <ConfirmClearDialog open={open} setOpen={setOpen} />
    </>
  )
})

interface ConfirmClearDialogProps {
  open: boolean
  setOpen: Dispatch<SetStateAction<boolean>>
}

const ConfirmClearDialog = memo(function ConfirmClearDialog({
  open,
  setOpen,
}: ConfirmClearDialogProps) {
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
      <DialogTitle>Confirm create new workflow</DialogTitle>
      <DialogContent>
        To create new workflow, current workflow will be cleared. The record
        will NOT be deleted if it has already been run. Are you sure?
      </DialogContent>
      <DialogActions>
        <Button variant="outlined" onClick={handleClose}>
          Cancel
        </Button>
        <Button variant="contained" onClick={handleClear}>
          OK
        </Button>
      </DialogActions>
    </Dialog>
  )
})

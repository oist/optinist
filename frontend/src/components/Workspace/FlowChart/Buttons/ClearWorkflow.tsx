import { Dispatch, memo, SetStateAction, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import DeleteIcon from "@mui/icons-material/Delete"
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
import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"

export const ClearWorkflowButton = memo(function ClearWorkflowButton() {
  const [open, setOpen] = useState(false)
  const openDialog = () => {
    setOpen(true)
  }

  const isPending = useSelector(selectPipelineIsStartedSuccess)

  return (
    <>
      <Tooltip title="Clear workflow">
        <IconButton
          disabled={!!isPending}
          onClick={openDialog}
          color={"primary"}
        >
          <DeleteIcon />
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
      <DialogTitle>Confirm clear workflow</DialogTitle>
      <DialogContent>
        Are you sure you want to clear the current workflow?
      </DialogContent>
      <DialogActions>
        <Button variant="outlined" onClick={handleClose}>
          Cancel
        </Button>
        <Button variant="contained" onClick={handleClear}>
          Clear
        </Button>
      </DialogActions>
    </Dialog>
  )
})

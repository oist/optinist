import { KeyboardEvent, memo } from "react"

import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Typography,
} from "@mui/material"

type ClearWorkflowIdDialogProps = {
  open: boolean
  handleOk: () => void
  handleCancel: () => void
}

export const ClearWorkflowIdDialog = memo(function ClearWorkflowIdDialog({
  open,
  handleOk,
  handleCancel,
}: ClearWorkflowIdDialogProps) {
  const handleClosePopup = (event: KeyboardEvent) => {
    if (event.key === "Escape") {
      handleCancel()
    }
  }
  return (
    <Dialog
      open={open}
      onClose={handleCancel}
      onKeyDown={handleClosePopup}
      fullWidth
    >
      <DialogTitle>Confirm running as new workflow</DialogTitle>
      <DialogContent>
        <Typography>
          {"With this operation, all algorithms should be run again."}
        </Typography>
        <Typography>
          {
            "So, this workflow will be run as new workflow with new ID (RUN ALL)."
          }
        </Typography>
        <Typography>{"Existing nodes are kept as they are."}</Typography>
        <Typography>{"Are you sure?"}</Typography>
      </DialogContent>
      <DialogActions>
        <Button onClick={handleCancel} variant="outlined">
          cancel
        </Button>
        <Button onClick={handleOk} variant="contained">
          ok
        </Button>
      </DialogActions>
    </Dialog>
  )
})

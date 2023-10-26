import { Dispatch, memo, SetStateAction, useContext, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar } from "notistack"

import ReplyIcon from "@mui/icons-material/Reply"
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Typography,
} from "@mui/material"
import IconButton from "@mui/material/IconButton"

import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable"
import { selectExperimentName } from "store/slice/Experiments/ExperimentsSelectors"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const ReproduceButton = memo(function ReproduceButton() {
  const [open, setOpen] = useState(false)
  const openDialog = () => {
    setOpen(true)
  }

  return (
    <>
      <IconButton onClick={openDialog}>
        <ReplyIcon />
      </IconButton>
      <ConfirmReprodceDialog open={open} setOpen={setOpen} />
    </>
  )
})

interface ConfirmReprodceDialogProps {
  open: boolean
  setOpen: Dispatch<SetStateAction<boolean>>
}

const ConfirmReprodceDialog = memo(function ConfirmReprodceDialog({
  open,
  setOpen,
}: ConfirmReprodceDialogProps) {
  const dispatch: AppDispatch = useDispatch()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = useContext(ExperimentUidContext)
  const workflowName = useSelector(selectExperimentName(uid))
  const { enqueueSnackbar } = useSnackbar()

  const handleClose = () => {
    setOpen(false)
  }

  const onClick = () => {
    if (workspaceId) {
      dispatch(reproduceWorkflow({ workspaceId, uid }))
        .unwrap()
        .then(() => {
          enqueueSnackbar("Successfully reproduced.", { variant: "success" })
          dispatch(reset())
        })
        .catch(() => {
          enqueueSnackbar("Failed to reproduce", { variant: "error" })
        })
    } else {
      enqueueSnackbar("Workspace id is missing", { variant: "error" })
    }
  }
  return (
    <Dialog open={open} onClose={handleClose}>
      <DialogTitle>Confirm reproduce workflow</DialogTitle>
      <DialogContent>
        <Typography>
          Reproduce <span style={{ fontWeight: "bold" }}>{workflowName}</span> (
          {uid})?
        </Typography>
      </DialogContent>
      <DialogActions>
        <Button variant="outlined" onClick={handleClose}>
          Cancel
        </Button>
        <Button variant="contained" onClick={onClick}>
          Reproduce
        </Button>
      </DialogActions>
    </Dialog>
  )
})

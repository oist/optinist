import { KeyboardEvent, memo, useContext, useState } from "react"
import { useSelector, useDispatch } from "react-redux"

import DeleteOutlineIcon from "@mui/icons-material/DeleteOutline"
import Button from "@mui/material/Button"
import Dialog from "@mui/material/Dialog"
import DialogActions from "@mui/material/DialogActions"
import DialogTitle from "@mui/material/DialogTitle"
import IconButton from "@mui/material/IconButton"

import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable"
import { deleteExperimentByUid } from "store/slice/Experiments/ExperimentsActions"
import { selectExperimentName } from "store/slice/Experiments/ExperimentsSelectors"
import {
  selectPipelineLatestUid,
  selectPipelineIsStartedSuccess,
} from "store/slice/Pipeline/PipelineSelectors"
import { clearCurrentPipeline } from "store/slice/Pipeline/PipelineSlice"
import { AppDispatch, RootState } from "store/store"

export const DeleteButton = memo(function DeleteButton() {
  const dispatch = useDispatch<AppDispatch>()
  const currentPipelineUid = useSelector(selectPipelineLatestUid)
  const uid = useContext(ExperimentUidContext)
  const isRunning = useSelector((state: RootState) => {
    const currentUid = selectPipelineLatestUid(state)
    const isPending = selectPipelineIsStartedSuccess(state)
    return uid === currentUid && isPending
  })
  const name = useSelector(selectExperimentName(uid))
  const [open, setOpen] = useState(false)

  const onClickOpen = () => {
    setOpen(true)
  }
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(deleteExperimentByUid(uid))
    uid === currentPipelineUid && dispatch(clearCurrentPipeline())
  }
  const handleClosePopup = (event: KeyboardEvent) => {
    if (event.key === "Escape") {
      onClickCancel()
    }
  }
  return (
    <>
      <IconButton onClick={onClickOpen} disabled={isRunning} color="error">
        <DeleteOutlineIcon />
      </IconButton>
      <Dialog open={open} onClose={onClickCancel} onKeyDown={handleClosePopup}>
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

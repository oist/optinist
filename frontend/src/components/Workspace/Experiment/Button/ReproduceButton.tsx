import { memo, useContext, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar } from "notistack"

import ReplyIcon from "@mui/icons-material/Reply"
import IconButton from "@mui/material/IconButton"

import { ConfirmDialog } from "components/common/ConfirmDialog"
import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable"
import { selectExperimentName } from "store/slice/Experiments/ExperimentsSelectors"
import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const ReproduceButton = memo(function ReproduceButton() {
  const [open, setOpen] = useState(false)
  const openDialog = () => {
    setOpen(true)
  }

  const isPending = useSelector(selectPipelineIsStartedSuccess)

  const dispatch: AppDispatch = useDispatch()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = useContext(ExperimentUidContext)
  const workflowName = useSelector(selectExperimentName(uid))
  const { enqueueSnackbar } = useSnackbar()
  const handleOk = () => {
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
    <>
      <IconButton onClick={openDialog} color="primary" disabled={!!isPending}>
        <ReplyIcon />
      </IconButton>
      <ConfirmDialog
        open={open}
        setOpen={setOpen}
        onConfirm={handleOk}
        title="Reproduce workflow?"
        content={`${workflowName} (${uid})`}
        confirmLabel="reproduce"
        iconType="info"
      />
    </>
  )
})

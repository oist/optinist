import { memo, useState } from "react"
import { useDispatch, useSelector } from "react-redux"

import { AddToPhotos } from "@mui/icons-material"
import { IconButton, Tooltip } from "@mui/material"

import { ConfirmDialog } from "components/common/ConfirmDialog"
import { clearFlowElements } from "store/slice/FlowElement/FlowElementSlice"
import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"

export const CreateWorkflowButton = memo(function CreateWorkflowButton() {
  const [open, setOpen] = useState(false)
  const dispatch = useDispatch()
  const isPending = useSelector(selectPipelineIsStartedSuccess)

  const openDialog = () => {
    setOpen(true)
  }
  const onConfirm = () => {
    dispatch(clearFlowElements())
  }

  return (
    <>
      <Tooltip title="Create new workflow">
        <IconButton onClick={openDialog} color="primary" disabled={!!isPending}>
          <AddToPhotos />
        </IconButton>
      </Tooltip>
      <ConfirmDialog
        open={open}
        setOpen={setOpen}
        onConfirm={onConfirm}
        title="Create new workflow?"
        content={`Current workflow will be cleared.
        If the workflow has already been run, the record will NOT be deleted.`}
        iconType="info"
      />
    </>
  )
})

import { memo } from "react"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar } from "notistack"

import AddchartIcon from "@mui/icons-material/Addchart"
import { IconButton, Tooltip } from "@mui/material"

import { getExperiments } from "store/slice/Experiments/ExperimentsActions"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { importSampleData } from "store/slice/Workflow/WorkflowActions"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

export const ImportSampleDataButton = memo(function ImportSampleDataButton() {
  const dispatch: AppDispatch = useDispatch()
  const { enqueueSnackbar } = useSnackbar()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const onClick = () => {
    if (typeof workspaceId === "number") {
      dispatch(importSampleData({ workspaceId }))
        .unwrap()
        .then(() => {
          enqueueSnackbar("Sample data import success", { variant: "success" })
          dispatch(reset())
          dispatch(getExperiments())
        })
        .catch(() => {
          enqueueSnackbar("Sample data import error", { variant: "error" })
        })
    }
  }
  return (
    <Tooltip title="Import sample data">
      <IconButton onClick={onClick}>
        <AddchartIcon style={{ color: "grey" }} />
      </IconButton>
    </Tooltip>
  )
})

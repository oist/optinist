import { memo } from "react"
import { useDispatch } from "react-redux"

import { useSnackbar } from "notistack"

import SdIcon from "@mui/icons-material/Sd"
import { IconButton, Tooltip } from "@mui/material"

import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { importSampleData } from "store/slice/Workflow/WorkflowActions"
import { AppDispatch } from "store/store"

export const ImportSampleDataButton = memo(function ImportSampleDataButton() {
  const dispatch: AppDispatch = useDispatch()
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    dispatch(importSampleData())
      .unwrap()
      .then(() => {
        enqueueSnackbar("Sample data import success", { variant: "success" })
        dispatch(reset())
      })
      .catch(() => {
        enqueueSnackbar("Sample data import error", { variant: "error" })
      })
  }

  return (
    <Tooltip title="Import sample data">
      <IconButton onClick={onClick}>
        <SdIcon color="primary" />
      </IconButton>
    </Tooltip>
  )
})

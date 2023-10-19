import React from "react"
import { useDispatch } from "react-redux"

import { useSnackbar } from "notistack"

import { UploadFile } from "@mui/icons-material"
import { IconButton, Tooltip } from "@mui/material"

import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { importWorkflowConfig } from "store/slice/Workflow/WorkflowActions"
import { AppDispatch } from "store/store"

export const ImportWorkflowConfigButton = React.memo(() => {
  const dispatch: AppDispatch = useDispatch()
  const inputRef = React.useRef<HTMLInputElement>(null)
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    inputRef.current?.click()
  }

  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    event.preventDefault()
    if (event.target.files != null && event.target.files[0] != null) {
      const file = event.target.files[0]
      const formData = new FormData()
      formData.append("file", file)
      dispatch(importWorkflowConfig({ formData }))
        .unwrap()
        .then(() => {
          enqueueSnackbar("Import success", { variant: "success" })
          dispatch(reset())
        })
        .catch(() => {
          enqueueSnackbar("Invalid yaml file", { variant: "error" })
        })
    }
  }

  return (
    <Tooltip title="Import workflow from config file">
      <IconButton onClick={onClick}>
        <UploadFile color="primary" />
        <input
          hidden
          ref={inputRef}
          type="file"
          accept=".yaml,.yml"
          onChange={onChange}
        />
      </IconButton>
    </Tooltip>
  )
})

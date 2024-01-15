import { memo, ChangeEvent, useRef } from "react"
import { useDispatch, useSelector } from "react-redux"

import { useSnackbar } from "notistack"

import { UploadFile } from "@mui/icons-material"
import { IconButton, Tooltip } from "@mui/material"

import { selectPipelineIsStartedSuccess } from "store/slice/Pipeline/PipelineSelectors"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { importWorkflowConfig } from "store/slice/Workflow/WorkflowActions"
import { AppDispatch } from "store/store"

export const ImportWorkflowConfigButton = memo(
  function ImportWorkflowConfigButton() {
    const dispatch: AppDispatch = useDispatch()
    const inputRef = useRef<HTMLInputElement>(null)
    const { enqueueSnackbar } = useSnackbar()

    const isPending = useSelector(selectPipelineIsStartedSuccess)

    const onClick = () => {
      inputRef.current?.click()
    }

    const onChange = (event: ChangeEvent<HTMLInputElement>) => {
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
      <Tooltip title="Import workflow yaml file">
        <IconButton onClick={onClick} color="primary" disabled={!!isPending}>
          <UploadFile />
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
  },
)

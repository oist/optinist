import { UploadFile } from '@mui/icons-material'
import { Button } from '@mui/material'
import { useSnackbar } from 'notistack'
import React from 'react'
import { useDispatch } from 'react-redux'
import { reset } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { importWorkflowConfig } from 'store/slice/Workflow/WorkflowActions'
import { AppDispatch } from 'store/store'

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
      formData.append('file', file)
      dispatch(importWorkflowConfig({ formData }))
        .unwrap()
        .then(() => {
          enqueueSnackbar('Import success', { variant: 'success' })
          dispatch(reset())
        })
        .catch(() => {
          enqueueSnackbar('Invalid yaml file', { variant: 'error' })
        })
    }
  }

  return (
    <Button
      variant="outlined"
      onClick={onClick}
      sx={{ margin: (theme) => theme.spacing(1) }}
      endIcon={<UploadFile />}
    >
      Import
      <input
        hidden
        ref={inputRef}
        type="file"
        accept=".yaml,.yml"
        onChange={onChange}
      />
    </Button>
  )
})

import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import IconButton from '@mui/material/IconButton'
import { useSnackbar } from 'notistack'
import { reproduceWorkflow } from 'store/slice/Workflow/WorkflowActions'
import { AppDispatch } from 'store/store'
import { ExperimentUidContext } from '../ExperimentTable'
import ReplyIcon from '@mui/icons-material/Reply'
import { reset } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { selectCurrentWorkspaceId } from 'store/slice/Workspace/WorkspaceSelector'

export const ReproduceButton = React.memo(() => {
  const dispatch: AppDispatch = useDispatch()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = React.useContext(ExperimentUidContext)
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    if (workspaceId) {
      dispatch(reproduceWorkflow({workspaceId, uid}))
      .unwrap()
      .then(() => {
        enqueueSnackbar('Successfully reproduced.', { variant: 'success' })
        dispatch(reset())
      })
      .catch(() => {
        enqueueSnackbar('Failed to reproduce', { variant: 'error' })
      })
    } else {
      enqueueSnackbar('Workspace id is missing', { variant: 'error' })
    }
  }
  return (
    <IconButton onClick={onClick}>
      <ReplyIcon />
    </IconButton>
  )
})

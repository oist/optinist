import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import IconButton from '@mui/material/IconButton'
import { useSnackbar } from 'notistack'
import { importExperimentByUid } from 'store/slice/Experiments/ExperimentsActions'
import { AppDispatch } from 'store/store'
import { ExperimentUidContext } from '../ExperimentTable'
import ReplyIcon from '@mui/icons-material/Reply'
import { reset } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { selectCurrentWorkspaceId } from 'store/slice/Workspace/WorkspaceSelector'

export const ImportButton = React.memo(() => {
  const dispatch: AppDispatch = useDispatch()
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const uid = React.useContext(ExperimentUidContext)
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    if (workspaceId) {
      dispatch(importExperimentByUid({workspaceId, uid}))
      .unwrap()
      .then(() => {
        enqueueSnackbar('Successfully imported.', { variant: 'success' })
        dispatch(reset())
      })
    } else {
      throw new Error('Workspace ID is missing.')
    }
  }
  return (
    <IconButton onClick={onClick}>
      <ReplyIcon />
    </IconButton>
  )
})

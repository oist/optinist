import React from 'react'
import { useDispatch } from 'react-redux'
import IconButton from '@mui/material/IconButton'
import { useSnackbar } from 'notistack'
import { importExperimentByUid } from 'store/slice/Experiments/ExperimentsActions'
import { AppDispatch } from 'store/store'
import { setRunBtnOption } from 'store/slice/Pipeline/PipelineSlice'
import { RUN_BTN_OPTIONS } from 'store/slice/Pipeline/PipelineType'
import { ExperimentUidContext } from '../ExperimentTable'
import ReplyIcon from '@mui/icons-material/Reply'

export const ImportButton = React.memo(() => {
  const dispatch: AppDispatch = useDispatch()
  const uid = React.useContext(ExperimentUidContext)
  const { enqueueSnackbar } = useSnackbar()

  const onClick = () => {
    dispatch(importExperimentByUid(uid))
      .unwrap()
      .then(() =>
        enqueueSnackbar('Successfully imported.', { variant: 'success' }),
      )
    dispatch(setRunBtnOption({ runBtnOption: RUN_BTN_OPTIONS.RUN_ALREADY }))
  }
  return (
    <IconButton onClick={onClick}>
      <ReplyIcon />
    </IconButton>
  )
})

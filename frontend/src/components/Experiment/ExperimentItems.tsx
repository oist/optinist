import React, { useState, useRef } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@mui/material/Button'
import IconButton from '@mui/material/IconButton'
import Dialog from '@mui/material/Dialog'
import DialogActions from '@mui/material/DialogActions'
import DialogTitle from '@mui/material/DialogTitle'
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline'
import GetAppIcon from '@mui/icons-material/GetApp'
import SimCardDownloadIcon from '@mui/icons-material/SimCardDownload'
import { useSnackbar } from 'notistack'
import { selectExperimentName } from 'store/slice/Experiments/ExperimentsSelectors'
import {
  deleteExperimentByUid,
  importExperimentByUid,
} from 'store/slice/Experiments/ExperimentsActions'
import { AppDispatch } from 'store/store'
import { setRunBtnOption } from 'store/slice/Pipeline/PipelineSlice'
import { RUN_BTN_OPTIONS } from 'store/slice/Pipeline/PipelineType'
import { ExperimentUidContext } from './ExperimentTable'
import { BASE_URL } from 'const/API'
import axios from 'axios'

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
      <GetAppIcon color="primary" />
    </IconButton>
  )
})

export const DeleteButton = React.memo(() => {
  const dispatch = useDispatch()
  const uid = React.useContext(ExperimentUidContext)

  const name = useSelector(selectExperimentName(uid))
  const [open, setOpen] = React.useState(false)

  const onClickOpen = () => {
    setOpen(true)
  }
  const onClickCancel = () => {
    setOpen(false)
  }
  const onClickOk = () => {
    setOpen(false)
    dispatch(deleteExperimentByUid(uid))
  }
  return (
    <>
      <IconButton onClick={onClickOpen}>
        <DeleteOutlineIcon color="error" />
      </IconButton>
      <Dialog open={open}>
        <DialogTitle>Are you sure you want to delete {name}?</DialogTitle>
        <DialogActions>
          <Button onClick={onClickCancel} variant="outlined" color="inherit">
            Cancel
          </Button>
          <Button onClick={onClickOk} variant="outlined" autoFocus>
            OK
          </Button>
        </DialogActions>
      </Dialog>
    </>
  )
})

export const DownloadButton = React.memo(() => {
  const uid = React.useContext(ExperimentUidContext)
  const ref = useRef<HTMLAnchorElement | null>(null)
  const [url, setFileUrl] = useState<string>()

  const onClick = async () => {
    try {
      const response = await axios.get(
        `${BASE_URL}/experiments/download/${uid}`,
        {
          responseType: 'blob',
        },
      )
      const url = URL.createObjectURL(new Blob([response.data]))
      setFileUrl(url)
      ref.current?.click()
      URL.revokeObjectURL(url)
    } catch (error) {
      throw new Error('Download Error')
    }
  }
  console.log(ref)

  return (
    <>
      <IconButton onClick={onClick}>
        <SimCardDownloadIcon color="primary" />
      </IconButton>
      <a href={url} download={`${uid}.nwb`} className="hidden" ref={ref} />
    </>
  )
})

import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import {
  selectPipelineIsCanceled,
  selectPipelineIsStartedSuccess,
  selectPipelineLatestUid,
  selectPipelineStatus,
  selectRunPostData,
} from './PipelineSelectors'
import { run, pollRunResult, runByCurrentUid } from './PipelineActions'
import { cancelPipeline } from './PipelineSlice'
import { AppDispatch } from 'store/store'
import { selectFilePathIsUndefined } from '../InputNode/InputNodeSelectors'
import { useSnackbar } from 'notistack'
import { RUN_STATUS } from './PipelineType'

const POLLING_INTERVAL = 5000

export type UseRunPipelineReturnType = ReturnType<typeof useRunPipeline>

export function useRunPipeline() {
  const dispatch: AppDispatch = useDispatch()
  const uid = useSelector(selectPipelineLatestUid)
  const isCanceled = useSelector(selectPipelineIsCanceled)
  const isStartedSuccess = useSelector(selectPipelineIsStartedSuccess)
  const filePathIsUndefined = useSelector(selectFilePathIsUndefined)
  const runPostData = useSelector(selectRunPostData)
  const handleRunPipeline = React.useCallback(() => {
    dispatch(run({ runPostData }))
  }, [dispatch, runPostData])
  const handleRunPipelineByUid = React.useCallback(() => {
    dispatch(runByCurrentUid({ runPostData }))
  }, [dispatch, runPostData])
  const handleCancelPipeline = React.useCallback(() => {
    if (uid != null) {
      dispatch(cancelPipeline())
    }
  }, [dispatch, uid])
  React.useEffect(() => {
    const intervalId = setInterval(() => {
      if (isStartedSuccess && !isCanceled && uid != null) {
        dispatch(pollRunResult({ uid: uid }))
      }
    }, POLLING_INTERVAL)
    return () => {
      clearInterval(intervalId)
    }
  }, [dispatch, uid, isCanceled, isStartedSuccess])
  const status = useSelector(selectPipelineStatus)
  const { enqueueSnackbar } = useSnackbar()
  // タブ移動による再レンダリングするたびにスナックバーが実行されてしまう挙動を回避するために前回の値を保持
  const [prevStatus, setPrevStatus] = React.useState(status)
  React.useEffect(() => {
    if (prevStatus !== status) {
      if (status === RUN_STATUS.FINISHED) {
        enqueueSnackbar('Finished', { variant: 'success' })
      } else if (status === RUN_STATUS.ABORTED) {
        enqueueSnackbar('Aborted', { variant: 'error' })
      } else if (status === RUN_STATUS.CANCELED) {
        enqueueSnackbar('Canceled', { variant: 'info' })
      }
      setPrevStatus(status)
    }
  }, [status, prevStatus, enqueueSnackbar])
  return {
    filePathIsUndefined,
    uid,
    status,
    handleRunPipeline,
    handleRunPipelineByUid,
    handleCancelPipeline,
  }
}

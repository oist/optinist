import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { selectRunPostData } from 'store/selectors/run/RunSelectors'
import {
  selectPipelineIsCanceled,
  selectPipelineIsStartedSuccess,
  selectPipelineLatestUid,
  selectPipelineRunBtn,
  selectPipelineStatus,
} from './PipelineSelectors'
import { run, pollRunResult, runByCurrentUid } from './PipelineActions'
import { cancelPipeline } from './PipelineSlice'
import { selectFilePathIsUndefined } from '../InputNode/InputNodeSelectors'
import { selectAlgorithmNodeNotExist } from '../AlgorithmNode/AlgorithmNodeSelectors'
import { useSnackbar } from 'notistack'
import { RUN_BTN_OPTIONS, RUN_STATUS } from './PipelineType'

const POLLING_INTERVAL = 5000

export type UseRunPipelineReturnType = ReturnType<typeof useRunPipeline>

export function useRunPipeline() {
  const dispatch = useDispatch()
  const uid = useSelector(selectPipelineLatestUid)
  const isCanceled = useSelector(selectPipelineIsCanceled)
  const isStartedSuccess = useSelector(selectPipelineIsStartedSuccess)
  const filePathIsUndefined = useSelector(selectFilePathIsUndefined)
  const algorithmNodeNotExist = useSelector(selectAlgorithmNodeNotExist)
  const runPostData = useSelector(selectRunPostData)
  const runButton = useSelector(selectPipelineRunBtn)
  const handleRunPipeline = React.useCallback(
    (name: string) => {
      dispatch(run({ runPostData: { name, ...runPostData, forceRunList: [] } }))
    },
    [dispatch, runPostData],
  )
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
        dispatch(
          pollRunResult({
            uid: uid,
            isReproduce: runButton === RUN_BTN_OPTIONS.RUN_ALREADY,
          }),
        )
      }
    }, POLLING_INTERVAL)
    return () => {
      clearInterval(intervalId)
    }
  }, [dispatch, uid, isCanceled, isStartedSuccess, runButton])
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
    algorithmNodeNotExist,
    uid,
    status,
    isStartedSuccess,
    handleRunPipeline,
    handleRunPipelineByUid,
    handleCancelPipeline,
  }
}

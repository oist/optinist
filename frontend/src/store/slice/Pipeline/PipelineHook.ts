import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import {
  selectPipelineIsCanceled,
  selectPipelineIsStartedSuccess,
  selectPipelineStatus,
  selectRunPostData,
} from './PipelineSelectors'
import { run, pollRunResult } from './PipelineActions'
import { cancelPipeline } from './PipelineSlice'
import { AppDispatch, RootState } from 'store/store'

const POLLING_INTERVAL = 5000

export function useRunPipeline() {
  const dispatch: AppDispatch = useDispatch()
  const [uid, setUid] = React.useState<null | string>()
  const isCanceled = useSelector((state: RootState) =>
    uid != null ? selectPipelineIsCanceled(uid)(state) : false,
  )
  const isStartedSuccess = useSelector((state: RootState) =>
    uid != null ? selectPipelineIsStartedSuccess(uid)(state) : false,
  )
  const runPostData = useSelector(selectRunPostData)
  const handleRunPipeline = React.useCallback(
    (newUid?: string) => {
      dispatch(run({ uid: newUid, runPostData }))
        .unwrap()
        .then((result) => setUid(result))
    },
    [dispatch, runPostData],
  )
  const handleCancelPipeline = React.useCallback(() => {
    if (uid != null) {
      dispatch(cancelPipeline(uid))
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
  const status = useSelector((state: RootState) =>
    uid != null ? selectPipelineStatus(uid)(state) : false,
  )
  return { uid, status, handleRunPipeline, handleCancelPipeline }
}

import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { selectRunPostData } from 'store/selectors/run/RunSelectors'
import {
  selectPipelineIsCanceled,
  selectPipelineIsStartedSuccess,
  selectPipelineLatestUid,
  selectPipelineStatus,
} from './PipelineSelectors'
import { AppDispatch } from 'store/store'
import { run, pollRunResult, runByCurrentUid } from './PipelineActions'
import { cancelPipeline } from './PipelineSlice'
import { selectFilePathIsUndefined } from '../InputNode/InputNodeSelectors'
import { selectAlgorithmNodeNotExist } from '../AlgorithmNode/AlgorithmNodeSelectors'
import { useSnackbar } from 'notistack'
import { RUN_STATUS } from './PipelineType'
import {
  getLastExperimentUid,
  importExperimentByUid,
} from '../Experiments/ExperimentsActions'

const POLLING_INTERVAL = 5000

export type UseRunPipelineReturnType = ReturnType<typeof useRunPipeline>

export function useRunPipeline() {
  const dispatch = useDispatch()
  const appDispatch: AppDispatch = useDispatch()
  const uid = useSelector(selectPipelineLatestUid)
  const isCanceled = useSelector(selectPipelineIsCanceled)
  const isStartedSuccess = useSelector(selectPipelineIsStartedSuccess)
  const filePathIsUndefined = useSelector(selectFilePathIsUndefined)
  const algorithmNodeNotExist = useSelector(selectAlgorithmNodeNotExist)
  const runPostData = useSelector(selectRunPostData)
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
    appDispatch(getLastExperimentUid())
      .unwrap()
      .then((uid) => {
        if (uid) {
          appDispatch(importExperimentByUid(uid))
            .unwrap()
            .then((_runPostData) => {
              dispatch(
                runByCurrentUid({
                  runPostData: { ...runPostData, ..._runPostData },
                }),
              )
            })
        }
      })
    //eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])
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
    algorithmNodeNotExist,
    uid,
    status,
    isStartedSuccess,
    handleRunPipeline,
    handleRunPipelineByUid,
    handleCancelPipeline,
  }
}

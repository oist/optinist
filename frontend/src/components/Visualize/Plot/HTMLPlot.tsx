import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { LinearProgress, Typography } from '@mui/material'

import { DisplayDataContext } from '../DataContext'
import {
  selectHTMLData,
  selectHTMLDataError,
  selectHTMLDataIsFulfilled,
  selectHTMLDataIsInitialized,
  selectHTMLDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getHTMLData } from 'store/slice/DisplayData/DisplayDataActions'

export const HTMLPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectHTMLDataIsPending(path))
  const isInitialized = useSelector(selectHTMLDataIsInitialized(path))
  const error = useSelector(selectHTMLDataError(path))
  const isFulfilled = useSelector(selectHTMLDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getHTMLData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <HTMLPlotImple />
  } else {
    return null
  }
})

const HTMLPlotImple = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const htmlData = useSelector(selectHTMLData(path))

  return <div dangerouslySetInnerHTML={{ __html: htmlData }} />
})

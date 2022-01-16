import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'

import { DisplayDataContext } from '../DisplayDataItem'
import {
  selectTimeSeriesData,
  selectTimeSeriesDataError,
  selectTimeSeriesDataIsFulfilled,
  selectTimeSeriesDataIsInitialized,
  selectTimeSeriesDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getTimeSeriesData } from 'store/slice/DisplayData/DisplayDataActions'
import { TimeSeriesData } from 'store/slice/DisplayData/DisplayDataType'
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'

export const TimeSeries = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectTimeSeriesDataIsPending(path))
  const isInitialized = useSelector(selectTimeSeriesDataIsInitialized(path))
  const error = useSelector(selectTimeSeriesDataError(path))
  const isFulfilled = useSelector(selectTimeSeriesDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getTimeSeriesData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <TimeSeriesImple />
  } else {
    return null
  }
})

const TimeSeriesImple = React.memo(() => {
  const { filePath: path, nodeId } = React.useContext(DisplayDataContext)
  const label = useSelector(selectNodeLabelById(nodeId))
  const timeSeriesData = useSelector(
    selectTimeSeriesData(path),
    timeSeriesDataEqualityFn,
  )
  const data = React.useMemo(() => {
    if (timeSeriesData == null) {
      return []
    }
    return Object.keys(timeSeriesData['0']).map((_, i) => {
      return {
        name: `(${i})`,
        x: Object.keys(timeSeriesData),
        y: Object.values(timeSeriesData).map((value) => value[i]),
        visible: i === 0 ? true : 'legendonly',
      }
    })
  }, [timeSeriesData])

  const layout = React.useMemo(
    () => ({
      title: label,
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      autosize: true,
      height: 300,
    }),
    [label],
  )

  const config = React.useMemo(
    () => ({
      displayModeBar: true,
    }),
    [],
  )
  return <PlotlyChart data={data} layout={layout} config={config} />
})

function timeSeriesDataEqualityFn(
  a: TimeSeriesData | undefined,
  b: TimeSeriesData | undefined,
) {
  if (a != null && b != null) {
    const aArray = Object.entries(a)
    const bArray = Object.entries(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every(([aKey, aValue], i) => {
          const [bKey, bValue] = bArray[i]
          return bKey === aKey && nestEqualityFun(bValue, aValue)
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}

function nestEqualityFun(
  a: {
    [key: number]: number
  },
  b: {
    [key: number]: number
  },
) {
  const aArray = Object.entries(a)
  const bArray = Object.entries(b)
  return (
    a === b ||
    (aArray.length === bArray.length &&
      aArray.every(([aKey, aValue], i) => {
        const [bKey, bValue] = bArray[i]
        return bKey === aKey && bValue === aValue
      }))
  )
}

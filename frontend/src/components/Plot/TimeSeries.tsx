import React from 'react'
import { useSelector, useDispatch } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import {
  algoNameByIdSelector,
  outputPathValueByIdSelector,
} from 'store/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import {
  timeSeriesDataSelector,
  timeSeriesDataIsLoadedSelector,
} from 'store/slice/PlotData/PlotDataSelector'
import { getTimeSeriesData } from 'store/slice/PlotData/PlotDataAction'
import { TimeSeriesData } from 'store/slice/PlotData/PlotDataType'

export const TimeSeries = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const dispatch = useDispatch()
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const isLoaded = useSelector(timeSeriesDataIsLoadedSelector(path ?? ''))
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getTimeSeriesData({ path }))
    }
  }, [dispatch, isLoaded, path])
  if (isLoaded) {
    return <TimeSeriesImple />
  } else {
    return null
  }
})

const TimeSeriesImple = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const name = useSelector(algoNameByIdSelector(nodeId))
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const timeSeriesData = useSelector(
    timeSeriesDataSelector(path ?? ''),
    timeSeriesDataEqualityFn,
  )
  const data = React.useMemo(() => {
    if (timeSeriesData == null) {
      return []
    }
    return Object.keys(timeSeriesData['0']).map((_, i) => {
      return {
        name: `${name}(${i})`,
        x: Object.keys(timeSeriesData),
        y: Object.values(timeSeriesData).map((value) => value[i]),
        visible: i == 0 ? true : 'legendonly',
      }
    })
  }, [timeSeriesData, name])

  const layout = React.useMemo(
    () => ({
      title: name,
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      autosize: true,
      height: 300,
    }),
    [name],
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

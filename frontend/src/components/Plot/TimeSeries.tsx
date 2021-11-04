import React from 'react'
import { useSelector, useDispatch } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import {
  algoNameByIdSelector,
  outputPathValueByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import {
  timeSeriesDataByKeySelector,
  timeSeriesDataIsLoadedByKeySelector,
} from 'redux/slice/PlotData/PlotDataSelector'
import { getTimeSeriesData } from 'redux/slice/PlotData/PlotDataAction'
import { toPlotDataKey } from 'redux/slice/PlotData/PlotDataUtils'

export const TimeSeries = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const dispatch = useDispatch()
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const isLoaded = useSelector(
    timeSeriesDataIsLoadedByKeySelector(toPlotDataKey(nodeId, outputKey)),
  )
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getTimeSeriesData({ nodeId, outputKey, path }))
    }
  }, [dispatch, isLoaded, nodeId, outputKey, path])
  if (isLoaded) {
    return <TimeSeriesImple />
  } else {
    return null
  }
})

const TimeSeriesImple = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const name = useSelector(algoNameByIdSelector(nodeId))
  const timeSeriesData = useSelector(
    timeSeriesDataByKeySelector(toPlotDataKey(nodeId, outputKey)),
  )
  if (timeSeriesData == null) {
    return null
  }
  const data = Object.keys(timeSeriesData.data['0']).map((_, i) => {
    return {
      name: `${name}(${i})`,
      x: Object.keys(timeSeriesData.data),
      y: Object.values(timeSeriesData.data).map((value) => value[i]),
    }
  })
  const layout = {
    title: name,
    margin: {
      t: 60, // top
      l: 50, // left
      b: 30, // bottom
    },
    autosize: true,
    height: 300,
  }
  const config = {
    displayModeBar: true,
  }
  return <PlotlyChart data={data} layout={layout} config={config} />
})

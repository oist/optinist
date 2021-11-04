import React from 'react'

import { useSelector, useDispatch } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import { OutputPlotContext } from 'App'
import {
  heatMapDataByKeySelector,
  heatMapDataIsLoadedByKeySelector,
} from 'redux/slice/PlotData/PlotDataSelector'
import { outputPathValueByIdSelector } from 'redux/slice/Algorithm/AlgorithmSelector'
import { getHeatMapData } from 'redux/slice/PlotData/PlotDataAction'
import { toPlotDataKey } from 'redux/slice/PlotData/PlotDataUtils'

export const HeatMap = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const dispatch = useDispatch()
  const path = useSelector(outputPathValueByIdSelector(nodeId, outputKey))
  const isLoaded = useSelector(
    heatMapDataIsLoadedByKeySelector(toPlotDataKey(nodeId, outputKey)),
  )
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getHeatMapData({ nodeId, outputKey, path }))
    }
  }, [dispatch, isLoaded, nodeId, outputKey, path])
  if (isLoaded) {
    return <HeatMapImple />
  } else {
    return null
  }
})
const HeatMapImple = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const currentOutputData = useSelector(
    heatMapDataByKeySelector(toPlotDataKey(nodeId, outputKey)),
  )
  if (currentOutputData == null) {
    return null
  }
  const data = [
    {
      z: currentOutputData.data,
      x: currentOutputData.data.map((_, i) => i),
      y: currentOutputData.data[0].map((_, i) => i),
      type: 'heatmap',
      hoverongaps: false,
    },
  ]

  const layout = {
    title: 'dummy',
    margin: {
      t: 60, // top
      l: 50, // left
      b: 30, // bottom
    },
    autosize: true,
    height: 350,
  }
  const config = {
    displayModeBar: true,
  }
  return <PlotlyChart data={data} layout={layout} config={config} />
})

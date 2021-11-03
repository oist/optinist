import React from 'react'
import { useSelector } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import {
  algoNameByIdSelector,
  outputDataByIdSelector,
  outputPathTypeByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import { toOutputDataId } from 'redux/slice/Algorithm/AlgorithmUtils'

export const TimeSeries = React.memo(() => {
  // const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  return <TimeSeriesImple />
})

const TimeSeriesImple = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const name = useSelector(algoNameByIdSelector(nodeId))
  const currentOutputData = useSelector(
    outputDataByIdSelector(toOutputDataId(nodeId, outputKey)),
  )
  if (currentOutputData == null) {
    return null
  }
  const data = Object.keys(currentOutputData.data['0']).map((_, i) => {
    return {
      name: `${name}(${i})`,
      x: Object.keys(currentOutputData.data),
      y: Object.values(currentOutputData.data).map((value) => value[i]),
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

import React from 'react'

import { useSelector } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import { outputDataByIdSelector } from 'redux/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import { toOutputDataId } from 'redux/slice/Algorithm/AlgorithmUtils'

export const HeatMap = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const currentOutputData = useSelector(
    outputDataByIdSelector(toOutputDataId(nodeId, outputKey)),
  )
  if (currentOutputData == null) {
    return null
  }

  const data = [
    {
      z: currentOutputData.data,
      x: Object.keys(currentOutputData.data).map((_, i) => i),
      y: Object.keys(currentOutputData.data[0]).map((_, i) => i),
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

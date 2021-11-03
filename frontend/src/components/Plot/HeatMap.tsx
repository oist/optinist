import React from 'react'

// import PlotlyChart from 'react-plotlyjs-ts'

// import React from 'react'
import { useSelector } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

import {
  algoNameByIdSelector,
  outputDataByIdSelector,
  outputPathTypeByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { OutputPlotContext } from 'App'
import { toOutputDataId } from 'redux/slice/Algorithm/AlgorithmUtils'

// export const TimeSeries = React.memo(() => {
//   // const { nodeId, outputKey } = React.useContext(OutputPlotContext)
//   return <TimeSeriesImple />
// })

export const HeatMap = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(OutputPlotContext)
  const currentOutputData = useSelector(
    outputDataByIdSelector(toOutputDataId(nodeId, outputKey)),
  )
  if (currentOutputData == null) {
    return null
  }
  console.log(currentOutputData.data)
  // var arr = Object.keys(currentOutputData.data).map((key1, _) => {
  //   Object.keys(currentOutputData.data).map((key2, _) => {
  //     return currentOutputData.data[key1][key2]
  //   })
  // })
  // console.log(arr)

  const data = [
    {
      z: currentOutputData.data,
      x: Object.keys(currentOutputData.data).map((_, i) => i),
      y: Object.keys(currentOutputData.data[0]).map((_, i) => i),
      type: 'heatmap',
      hoverongaps: false,
    },
  ]
  console.log(data)
  // const data = [
  //   {
  //     z: [
  //       [1, 20, 30, 50, 1],
  //       [20, 1, 60, 80, 30],
  //       [30, 60, 1, -10, 20],
  //     ],
  //     x: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
  //     y: ['Morning', 'Afternoon', 'Evening'],
  //     type: 'heatmap',
  //     hoverongaps: false,
  //   },
  // ]
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

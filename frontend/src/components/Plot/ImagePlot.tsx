import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'

export const ImagePlot = React.memo(() => {
  return <ImagePlotImple />
})

const ImagePlotImple = React.memo(() => {
  const data = [
    {
      z: [
        [0, 0, 0],
        [0, 10, 10, 0],
        [0, 100, 100, 0],
        [0, 0, 0, 0],
      ],
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

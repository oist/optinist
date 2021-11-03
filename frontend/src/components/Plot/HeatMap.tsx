import React from 'react'

import PlotlyChart from 'react-plotlyjs-ts'

export const HeatMap = React.memo(() => {
  const data = [
    {
      z: [
        [1, 20, 30, 50, 1],
        [20, 1, 60, 80, 30],
        [30, 60, 1, -10, 20],
      ],
      x: ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday'],
      y: ['Morning', 'Afternoon', 'Evening'],
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

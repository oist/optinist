import React from 'react'
import { Line } from 'react-chartjs-2'
import { useSelector, useDispatch } from 'react-redux'
import { getAlgoOutputData } from 'redux/slice/Algorithm/AlgorithmAction'
import {
  currentAlgoNameByIdSelector,
  currentOutputDataSelector,
  outputDataIsLoadedByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { NodeIdContext } from 'App'

export const PlotOutput = React.memo(function PlotOutput() {
  const nodeId = React.useContext(NodeIdContext)
  const dispatch = useDispatch()
  const name = useSelector(currentAlgoNameByIdSelector(nodeId))
  const isLoaded = useSelector(outputDataIsLoadedByIdSelector(nodeId))
  React.useEffect(() => {
    if (!isLoaded && name) {
      dispatch(getAlgoOutputData({ id: nodeId, name }))
    }
  }, [isLoaded, nodeId])
  if (isLoaded) {
    return (
      <div>
        {name}
        <Chart />
      </div>
    )
  } else {
    return null
  }
})

const Chart = React.memo(() => {
  const nodeId = React.useContext(NodeIdContext)
  const name = useSelector(currentAlgoNameByIdSelector(nodeId))
  const currentOutputData = useSelector(currentOutputDataSelector(nodeId))
  if (currentOutputData == null) {
    return null
  }
  const data = {
    labels: Object.keys(currentOutputData.data),
    // .filter((_, i) => i % 50 == 0), // chart.jsが重いので間引く
    datasets: Object.keys(currentOutputData.data['0']).map((_, i) => {
      const color = '#' + Math.floor(Math.random() * 16777215).toString(16)
      return {
        label: `${name}(${i})`,
        data: Object.values(currentOutputData.data)
          // .filter((_, i) => i % 50 == 0)
          .map((value) => value[i]),
        backgroundColor: color,
        borderColor: color,
        borderWidth: 1,
      }
    }),
  }
  const options = {
    radius: 2,
    scales: {
      y: {
        beginAtZero: true,
      },
    },
  }
  return <Line data={data} width={100} height={60} options={options} />
})

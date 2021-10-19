import React from 'react'
import { Line } from 'react-chartjs-2'
import { useSelector, useDispatch } from 'react-redux'
import { RootState } from 'redux/store'
import { getAlgoOutputData } from 'redux/slice/Algorithm/AlgorithmAction'
import {
  currentAlgoIdSelector,
  currentAlgoNameSelector,
  currentOutputDataSelector,
  outputDataIsLoadedByIdSelector,
  selectedOutputPathSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'

export const PlotOutput = React.memo(function PlotOutput() {
  const dispatch = useDispatch()
  const id = useSelector(currentAlgoIdSelector)
  const name = useSelector(currentAlgoNameSelector)
  const isPlotData = useSelector((state: RootState) => {
    const isImage = selectedOutputPathSelector(id)(state)?.isImage
    return isImage != null ? !isImage : false
  })
  const isLoaded = useSelector(outputDataIsLoadedByIdSelector(id))
  React.useEffect(() => {
    if (isPlotData && !isLoaded && name) {
      dispatch(getAlgoOutputData({ id, name }))
    }
  }, [isPlotData, isLoaded, id])
  if (isPlotData && isLoaded) {
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
  const name = useSelector(currentAlgoNameSelector)
  const currentOutputData = useSelector(currentOutputDataSelector)
  if (currentOutputData == null) {
    return null
  }
  const data = {
    labels: Object.keys(currentOutputData.data).filter((_, i) => i % 50 == 0), // chart.jsが重いので間引く
    datasets: Object.keys(currentOutputData.data['0']).map((_, i) => {
      const color = '#' + Math.floor(Math.random() * 16777215).toString(16)
      return {
        label: `${name}(${i})`,
        data: Object.values(currentOutputData.data)
          .filter((_, i) => i % 50 == 0)
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
      yAxes: [
        {
          ticks: {
            beginAtZero: true,
          },
        },
      ],
    },
  }
  return <Line data={data} width={100} height={60} options={options} />
})

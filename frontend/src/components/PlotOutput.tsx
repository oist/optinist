import React from 'react'
import { Bar } from 'react-chartjs-2'
import { useSelector, useDispatch } from 'react-redux'
import { getAlgoOutputData } from 'redux/slice/Algorithm/AlgorithmAction'
import {
  currentOutputDataSelector,
  currentAlgoIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { RootState } from 'redux/store'

export const PlotOutput = React.memo(function PlotOutput() {
  const dispatch = useDispatch()
  const currentOutputDataIsLoaded = useSelector(
    (state: RootState) => currentOutputDataSelector(state) != null,
  )
  const currentAlgoId = useSelector(currentAlgoIdSelector)
  React.useEffect(() => {
    if (!currentOutputDataIsLoaded) {
      dispatch(getAlgoOutputData({ id: currentAlgoId }))
    }
  }, [currentAlgoId, currentOutputDataIsLoaded, dispatch])
  if (currentOutputDataIsLoaded) {
    return <Chart />
  } else {
    return null
  }
})

const Chart = React.memo(() => {
  const currentAlgoId = useSelector(currentAlgoIdSelector)
  const currentOutputData = useSelector(currentOutputDataSelector)
  if (currentOutputData == null) {
    return null
  }
  const data = {
    labels: currentOutputData.data.map((data) => data.x),
    datasets: [
      {
        label: 'Dataset(' + currentAlgoId + ')',
        data: currentOutputData.data?.map((data) => data.y),
        backgroundColor: 'rgb(255, 99, 132, 0.5)',
        borderColor: 'rgb(255, 99, 132)',
        borderWidth: 1,
      },
    ],
  }

  const options = {
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
  return <Bar data={data} width={100} height={40} options={options} />
})

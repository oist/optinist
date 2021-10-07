import React from 'react'
import { Bar } from 'react-chartjs-2'
import { useSelector, useDispatch } from 'react-redux'
import { getOutputData } from 'redux/slice/Output/OutputAction'
import {
  currentOutputDataSelector,
  currentOutputIdSelector,
} from 'redux/slice/Output/OutputSelector'
import { RootState } from 'redux/store'

export const PlotOutput = React.memo(function PlotOutput() {
  const dispatch = useDispatch()
  const currentOutputDataIsLoaded = useSelector(
    (state: RootState) => currentOutputDataSelector(state) != null,
  )
  const currentOutputId = useSelector(currentOutputIdSelector)
  React.useEffect(() => {
    if (!currentOutputDataIsLoaded) {
      dispatch(getOutputData({ id: currentOutputId }))
    }
  }, [currentOutputId, currentOutputDataIsLoaded, dispatch])
  if (currentOutputDataIsLoaded) {
    return <Chart />
  } else {
    return null
  }
})

const Chart = React.memo(() => {
  const currentOutputId = useSelector(currentOutputIdSelector)
  const currentOutputData = useSelector(currentOutputDataSelector)
  if (currentOutputData == null) {
    return null
  }
  const data = {
    labels: currentOutputData.map((data) => data.x),
    datasets: [
      {
        label: 'Dataset(' + currentOutputId + ')',
        data: currentOutputData?.map((data) => data.y),
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

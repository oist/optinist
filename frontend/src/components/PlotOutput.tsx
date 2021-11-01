import React from 'react'
import { useSelector, useDispatch } from 'react-redux'

import PlotlyChart from 'react-plotlyjs-ts'

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
    return <Chart />
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

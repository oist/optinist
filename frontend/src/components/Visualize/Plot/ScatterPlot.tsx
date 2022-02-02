import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'

import { DisplayDataContext } from '../DataContext'
import {
  selectScatterData,
  selectScatterDataError,
  selectScatterDataIsFulfilled,
  selectScatterDataIsInitialized,
  selectScatterDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getScatterData } from 'store/slice/DisplayData/DisplayDataActions'
import { ScatterData } from 'store/slice/DisplayData/DisplayDataType'

export const ScatterPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const dispatch = useDispatch()
  const isPending = useSelector(selectScatterDataIsPending(path))
  const isInitialized = useSelector(selectScatterDataIsInitialized(path))
  const error = useSelector(selectScatterDataError(path))
  const isFulfilled = useSelector(selectScatterDataIsFulfilled(path))
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getScatterData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <ScatterPlotImple />
  } else {
    return null
  }
})

const ScatterPlotImple = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)

  const scatterData = useSelector(
    selectScatterData(path),
    scatterDataEqualityFn,
  )

  console.log(scatterData)
  const data = React.useMemo(
    () => [
      {
        x: Object.values(scatterData[0]),
        y: Object.values(scatterData[1]),
        type: 'scatter',
        mode: 'markers', //'markers+text',
        text: Object.keys(scatterData[0]),
        textposition: 'top center',
        textfont: {
          family: 'Raleway, sans-serif',
        },
        marker: {
          size: 5,
          color: Object.keys(scatterData[0]),
        },
      },
    ],
    [scatterData],
  )

  const layout = React.useMemo(
    () => ({
      title: path.split('/').reverse()[0],
      margin: {
        t: 60, // top
        l: 50, // left
        b: 30, // bottom
      },
      dragmode: 'pan',
      autosize: true,
    }),
    [],
  )

  const config = {
    displayModeBar: true,
    scrollZoom: true,
  }

  return <PlotlyChart data={data} layout={layout} config={config} />
})

function scatterDataEqualityFn(
  a: ScatterData | undefined,
  b: ScatterData | undefined,
) {
  if (a != null && b != null) {
    const aArray = Object.entries(a)
    const bArray = Object.entries(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every(([aKey, aValue], i) => {
          const [bKey, bValue] = bArray[i]
          return bKey === aKey // && nestEqualityFun(bValue, aValue)
        }))
    )
  } else {
    return a === undefined && b === undefined
  }
}

// function nestEqualityFun(
//   a: {
//     [key: number]: number
//   },
//   b: {
//     [key: number]: number
//   },
// ) {
//   const aArray = Object.entries(a)
//   const bArray = Object.entries(b)
//   return (
//     a === b ||
//     (aArray.length === bArray.length &&
//       aArray.every(([aKey, aValue], i) => {
//         const [bKey, bValue] = bArray[i]
//         return bKey === aKey && bValue === aValue
//       }))
//   )
// }

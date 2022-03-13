import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@mui/material'

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
import {
  selectScatterItemXIndex,
  selectScatterItemYIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

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
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  const scatterData = useSelector(
    selectScatterData(path),
    scatterDataEqualityFn,
  )

  const xIndex = useSelector(selectScatterItemXIndex(itemId))
  const yIndex = useSelector(selectScatterItemYIndex(itemId))
  const maxIndex = Object.keys(scatterData).length - 1

  const data = React.useMemo(
    () => [
      {
        x:
          xIndex < maxIndex
            ? Object.values(scatterData[xIndex])
            : Object.values(scatterData[maxIndex]),
        y:
          yIndex < maxIndex
            ? Object.values(scatterData[yIndex])
            : Object.values(scatterData[maxIndex]),
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
    [scatterData, xIndex, yIndex, maxIndex],
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
      xaxis: {
        title: {
          text: `x: ${
            xIndex < Object.keys(scatterData).length ? xIndex : maxIndex
          }`,
          font: {
            family: 'Courier New, monospace',
            size: 18,
            color: '#7f7f7f',
          },
        },
      },
      yaxis: {
        title: {
          text: `y: ${
            yIndex < Object.keys(scatterData).length ? yIndex : maxIndex
          }`,
          font: {
            family: 'Courier New, monospace',
            size: 18,
            color: '#7f7f7f',
          },
        },
      },
    }),
    [xIndex, yIndex, maxIndex, path, scatterData],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
  }

  const ref = React.useRef<HTMLDivElement>(null)
  const plotlyHeight = ref.current?.getBoundingClientRect().height

  useEffect(() => {
    const height =
      ref.current?.getElementsByClassName('main-svg')[0].clientHeight
    const plotContainer = (
      ref.current?.getElementsByClassName(
        'plot-container',
      ) as HTMLCollectionOf<HTMLElement>
    )[0]

    if (height !== undefined && plotContainer !== undefined) {
      plotContainer.style.height = `${height}px`
    }
  }, [plotlyHeight])

  return (
    <div ref={ref}>
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

function scatterDataEqualityFn(
  a: ScatterData | undefined,
  b: ScatterData | undefined,
) {
  if (a != null && b != null) {
    const aArray = Object.keys(a)
    const bArray = Object.keys(b)
    return (
      a === b ||
      (aArray.length === bArray.length &&
        aArray.every((aKey, i) => {
          const bKey = bArray[i]
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

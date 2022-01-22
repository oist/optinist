import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import { DisplayDataContext } from '../DataContext'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectRoiData,
  selectRoiDataError,
  selectRoiDataIsFulfilled,
  selectRoiDataIsInitialized,
  selectRoiDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { LinearProgress, Typography } from '@material-ui/core'
import { getRoiData } from 'store/slice/DisplayData/DisplayDataActions'
import { selectRoiItemColors } from 'store/slice/VisualizeItem/VisualizeItemSelectors'

export const RoiPlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const isPending = useSelector(selectRoiDataIsPending(path))
  const isInitialized = useSelector(selectRoiDataIsInitialized(path))
  const isFulfilled = useSelector(selectRoiDataIsFulfilled(path))
  const error = useSelector(selectRoiDataError(path))

  const dispatch = useDispatch()
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getRoiData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <RoiPlotImple />
  } else {
    return null
  }
})

const RoiPlotImple = React.memo<{}>(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const imageData = useSelector(selectRoiData(path), imageDataEqualtyFn)

  // const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  // const showline = useSelector(selectImageItemShowLine(itemId))
  // const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  // const showgrid = useSelector(selectImageItemShowGrid(itemId))
  // const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectRoiItemColors(itemId))

  const data = React.useMemo(
    () => [
      {
        z: imageData,
        type: 'heatmap',
        name: 'images',
        colorscale: colorscale.map((value) => {
          let offset: number = parseFloat(value.offset)
          const offsets: number[] = colorscale.map((v) => {
            return parseFloat(v.offset)
          })
          // plotlyは端[0.0, 1.0]がないとダメなので、その設定
          if (offset === Math.max(...offsets)) {
            offset = 1.0
          }
          if (offset === Math.min(...offsets)) {
            offset = 0.0
          }
          return [offset, value.rgb]
        }),
        hoverongaps: false,
        // showscale: showscale,
        // zsmooth: zsmooth, // ['best', 'fast', false]
        zsmooth: false,
        showlegend: true,
      },
    ],
    [imageData, colorscale],
  )

  const layout = React.useMemo(
    () => ({
      title: path.split('/').reverse()[0],
      margin: {
        t: 30, // top
        l: 120, // left
        b: 30, // bottom
      },
      dragmode: 'pan',
      xaxis: {
        autorange: true,
        // showgrid: showgrid,
        // showline: showline,
        zeroline: false,
        autotick: true,
        ticks: '',
        // showticklabels: showticklabels,
      },
      yaxis: {
        autorange: 'reversed',
        // showgrid: showgrid,
        // showline: showline,
        zeroline: false,
        autotick: true, // todo
        ticks: '',
        // showticklabels: showticklabels, // todo
      },
    }),
    [path],
  )
  const config = {
    displayModeBar: true,
    scrollZoom: true,
  }
  return (
    <div className="imagePlotChart">
      <PlotlyChart data={data} layout={layout} config={config} />
    </div>
  )
})

function imageDataEqualtyFn(
  a: number[][] | undefined,
  b: number[][] | undefined,
) {
  if (a != null && b != null) {
    return twoDimarrayEqualityFn(a, b)
  } else {
    return a === undefined && b === undefined
  }
}

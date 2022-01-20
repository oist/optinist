import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import { DisplayDataContext } from '../DisplayDataItem'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectRoiActivateData,
  selectRoiData,
  selectRoiDataError,
  selectRoiDataIsFulfilled,
  selectRoiDataIsInitialized,
  selectRoiDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import {
  Button,
  LinearProgress,
  MobileStepper,
  Typography,
  useTheme,
} from '@material-ui/core'
import { getRoiData } from 'store/slice/DisplayData/DisplayDataActions'

export const RoiPlot = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
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
  const imageData = useSelector(selectRoiActivateData(path), imageDataEqualtyFn)
  console.log(imageData)

  // const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  // const showline = useSelector(selectImageItemShowLine(itemId))
  // const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  // const showgrid = useSelector(selectImageItemShowGrid(itemId))
  // const showscale = useSelector(selectImageItemShowScale(itemId))
  // const colorscale = useSelector(selectImageItemColors(itemId))

  const data = React.useMemo(
    () => [
      {
        z: imageData,
        type: 'heatmap',
        name: 'images',
        colorscale: [
          [0, 'rgb(255, 255, 255)'],
          [1.0, 'rgb(0, 0, 0)'],
        ],
        hoverongaps: false,
        // showscale: showscale,
        // zsmooth: zsmooth, // ['best', 'fast', false]
        showlegend: true,
      },
    ],
    [imageData],
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
    [],
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

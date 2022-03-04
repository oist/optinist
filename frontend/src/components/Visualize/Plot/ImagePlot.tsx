import React, { useEffect } from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import {
  Button,
  LinearProgress,
  MobileStepper,
  Typography,
  useTheme,
} from '@mui/material'
import KeyboardArrowLeft from '@mui/icons-material/KeyboardArrowLeft'
import KeyboardArrowRight from '@mui/icons-material/KeyboardArrowRight'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import { DisplayDataContext } from '../DataContext'

import {
  selectImageDataError,
  selectImageDataIsInitialized,
  selectImageDataIsPending,
  selectImageDataIsFulfilled,
  selectActiveImageData,
  selectRoiData,
  selectImageDataEndIndex,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import {
  getImageData,
  getRoiData,
  getTimeSeriesDataById,
} from 'store/slice/DisplayData/DisplayDataActions'
import {
  selectImageItemShowticklabels,
  selectImageItemZsmooth,
  selectImageItemShowLine,
  selectImageItemShowGrid,
  selectImageItemShowScale,
  selectImageItemColors,
  selectImageItemActiveIndex,
  selectImageItemStartIndex,
  selectImageItemEndIndex,
  selectRoiItemFilePath,
  selectMultiPlotTimeSeriesItemFilepath,
  selectMultiPlotTimeSeriesItemDisplayNumbers,
  selectRoiItemIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  decrementImageActiveIndex,
  incrementImageActiveIndex,
  setTimeSeriesItemDisplayNumbers,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { Datum, LayoutAxis, PlotData } from 'plotly.js'
import createColormap from 'colormap'

export const ImagePlot = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)

  const startIndex = useSelector(selectImageItemStartIndex(itemId))
  const endIndex = useSelector(selectImageItemEndIndex(itemId))
  const isPending = useSelector(selectImageDataIsPending(path))
  const isInitialized = useSelector(selectImageDataIsInitialized(path))
  const isFulfilled = useSelector(selectImageDataIsFulfilled(path))
  const error = useSelector(selectImageDataError(path))

  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))

  const dispatch = useDispatch()
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(
        getImageData({
          path,
          startIndex: startIndex ?? 1,
          endIndex: endIndex ?? 10,
        }),
      )
    }
    if (roiFilePath != null) {
      dispatch(getRoiData({ path: roiFilePath }))
    }
  }, [dispatch, isInitialized, path, startIndex, endIndex, roiFilePath])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <ImagePlotImple />
  } else {
    return null
  }
})

const ImagePlotImple = React.memo(() => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const itemStartIndex = useSelector(selectImageItemStartIndex(itemId))
  const itemEndIndex = useSelector(selectImageItemEndIndex(itemId))
  const endIndex = useSelector(selectImageDataEndIndex(path))
  const activeIndex = useSelector(selectImageItemActiveIndex(itemId))
  const dispatch = useDispatch()
  const handleNext = () => dispatch(incrementImageActiveIndex({ itemId }))
  const handleBack = () => dispatch(decrementImageActiveIndex({ itemId }))
  const theme = useTheme()
  return (
    <>
      <MobileStepper
        steps={itemEndIndex}
        position="static"
        variant="text"
        activeStep={activeIndex + itemStartIndex - 1}
        nextButton={
          <Button
            size="small"
            onClick={handleNext}
            disabled={activeIndex === (endIndex ?? 0)}
          >
            <Typography>Next</Typography>
            {theme.direction === 'rtl' ? (
              <KeyboardArrowLeft />
            ) : (
              <KeyboardArrowRight />
            )}
          </Button>
        }
        backButton={
          <Button
            size="small"
            onClick={handleBack}
            disabled={activeIndex === 0}
          >
            {theme.direction === 'rtl' ? (
              <KeyboardArrowRight />
            ) : (
              <KeyboardArrowLeft />
            )}
            <Typography>Back</Typography>
          </Button>
        }
      />
      <ImagePlotChart activeIndex={activeIndex} />
    </>
  )
})

const ImagePlotChart = React.memo<{
  activeIndex: number
}>(({ activeIndex }) => {
  const { filePath: path, itemId } = React.useContext(DisplayDataContext)
  const imageData = useSelector(
    selectActiveImageData(path, activeIndex),
    imageDataEqualtyFn,
  )
  const roiFilePath = useSelector(selectRoiItemFilePath(itemId))
  const roiData = useSelector(
    (state: RootState) =>
      roiFilePath != null ? selectRoiData(roiFilePath)(state) : [],
    imageDataEqualtyFn,
  )

  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))

  const timeDataMaxIndex = useSelector(selectRoiItemIndex(itemId, roiFilePath))

  const colorscaleRoi = createColormap({
    colormap: 'jet',
    nshades: 100, //timeDataMaxIndex >= 6 ? timeDataMaxIndex : 6,
    format: 'hex',
    alpha: 1,
  })

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
        showscale: showscale,
        zsmooth: zsmooth, // ['best', 'fast', false]
        // showlegend: true,
      },
      {
        z: roiData,
        type: 'heatmap',
        name: 'roi',
        hovertemplate: 'cell id: %{z}',
        colorscale: [...Array(timeDataMaxIndex)].map((_, i) => {
          const new_i = Math.floor((i % 10) * 10 + i / 10)
          const offset = i / (timeDataMaxIndex - 1)
          const rgb = colorscaleRoi[new_i]
          return [offset, rgb]
        }),
        zmin: 1,
        zmax: timeDataMaxIndex,
        hoverongaps: false,
        zsmooth: false,
        showscale: false,
      },
    ],
    [
      imageData,
      roiData,
      zsmooth,
      showscale,
      colorscale,
      colorscaleRoi,
      timeDataMaxIndex,
    ],
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
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true,
        ticks: '',
        showticklabels: showticklabels,
      },
      yaxis: {
        automargin: true,
        autorange: 'reversed',
        showgrid: showgrid,
        showline: showline,
        zeroline: false,
        autotick: true, // todo
        ticks: '',
        showticklabels: showticklabels, // todo
      },
      // annotations: [
      //   {
      //     xref: 'x1',
      //     yref: 'y1',
      //     x: 10,
      //     y: 10,
      //     text: "aaa",
      //     showarrow: false,
      //     font: {
      //       family: 'Arial',
      //       size: 30,
      //       color: 'rgb(255, 255, 255)'
      //     }
      //   }
      // ],
    }),
    [path, showgrid, showline, showticklabels],
  )

  const config = {
    displayModeBar: true,
    responsive: true,
    height: '100%',
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
  }, [plotlyHeight, activeIndex])

  const dispatch = useDispatch()
  const timeSeriesFilePath = useSelector(
    selectMultiPlotTimeSeriesItemFilepath(itemId),
  )
  const displayNumbers = useSelector(
    selectMultiPlotTimeSeriesItemDisplayNumbers(itemId),
  )

  const onClick = (event: any) => {
    const points: PlotDatum = event.points[0]
    if (
      timeSeriesFilePath !== null &&
      displayNumbers !== null &&
      points.curveNumber >= 1 &&
      !displayNumbers.includes(points.z - 1)
    ) {
      const newValue = [...displayNumbers, points.z - 1]
      dispatch(
        setTimeSeriesItemDisplayNumbers({
          itemId,
          displayNumbers: newValue,
        }),
      )
      dispatch(
        getTimeSeriesDataById({
          path: timeSeriesFilePath,
          index: points.z - 1,
        }),
      )
    }
  }

  return (
    <div ref={ref}>
      <PlotlyChart
        data={data}
        layout={layout}
        config={config}
        onClick={onClick}
      />
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

interface PlotDatum {
  curveNumber: number
  data: PlotData
  customdata: Datum
  pointIndex: number
  pointNumber: number
  x: Datum
  xaxis: LayoutAxis
  y: Datum
  yaxis: LayoutAxis
  z: number
}

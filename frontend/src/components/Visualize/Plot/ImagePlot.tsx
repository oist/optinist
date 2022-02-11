import React, { useEffect } from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import {
  Button,
  LinearProgress,
  MobileStepper,
  Typography,
  useTheme,
} from '@material-ui/core'
import KeyboardArrowLeft from '@material-ui/icons/KeyboardArrowLeft'
import KeyboardArrowRight from '@material-ui/icons/KeyboardArrowRight'

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
  selectRoiItemColors,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  decrementImageActiveIndex,
  incrementImageActiveIndex,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'

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
  const colorscaleRoi = useSelector(selectRoiItemColors(itemId))

  const showticklabels = useSelector(selectImageItemShowticklabels(itemId))
  const showline = useSelector(selectImageItemShowLine(itemId))
  const zsmooth = useSelector(selectImageItemZsmooth(itemId))
  const showgrid = useSelector(selectImageItemShowGrid(itemId))
  const showscale = useSelector(selectImageItemShowScale(itemId))
  const colorscale = useSelector(selectImageItemColors(itemId))

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
        colorscale: colorscaleRoi.map((value) => {
          let offset: number = parseFloat(value.offset)
          const offsets: number[] = colorscaleRoi.map((v) => {
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
        zsmooth: false, // ['best', 'fast', false]
        showlegend: true,
        showscale: false,
      },
    ],
    [imageData, roiData, zsmooth, showscale, colorscale, colorscaleRoi],
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
    }),
    [path, showgrid, showline, showticklabels],
  )
  const config = {
    displayModeBar: true,
    // scrollZoom: true,
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

  return (
    <div ref={ref}>
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

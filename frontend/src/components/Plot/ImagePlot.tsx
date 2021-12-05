import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { useSelector, useDispatch } from 'react-redux'
import { ImageDataContext } from 'App'
import { outputPathValueByIdSelector } from 'store/slice/Algorithm/AlgorithmSelector'
import {
  imageDataActiveIndexselector,
  activeImageDataSelector,
  imageDataIsLoadedSelector,
  imageDataMaxIndexSelector,
} from 'store/slice/PlotData/PlotDataSelector'
import { getImageData } from 'store/slice/PlotData/PlotDataAction'
import { RootState } from 'store/store'
import {
  imageIsUploadedByIdSelector,
  imageTiffPathByIdSelector,
  imageMaxIndexByIdSelector,
} from 'store/slice/ImageFile/ImageFileSelector'
import {
  Button,
  LinearProgress,
  MobileStepper,
  Typography,
  useTheme,
} from '@material-ui/core'
import KeyboardArrowLeft from '@material-ui/icons/KeyboardArrowLeft'
import KeyboardArrowRight from '@material-ui/icons/KeyboardArrowRight'
import {
  decrementImageActiveIndex,
  incrementImageActiveIndex,
} from 'store/slice/PlotData/PlotData'
import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

export const ImagePlot = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(ImageDataContext)
  const isUploaded = useSelector(imageIsUploadedByIdSelector(nodeId))
  if (outputKey != null || isUploaded === true) {
    return <ImagePlotContainer />
  } else {
    return null
  }
})

const ImagePlotContainer = React.memo(() => {
  const { nodeId, outputKey } = React.useContext(ImageDataContext)
  const dispatch = useDispatch()
  const path = useSelector((state: RootState) => {
    if (outputKey != null) {
      // Algoの出力データの場合
      return outputPathValueByIdSelector(nodeId, outputKey)(state)
    } else {
      // 画像ファイルの入力データの場合
      return imageTiffPathByIdSelector(nodeId)(state)
    }
  })
  const maxIndex = useSelector(imageMaxIndexByIdSelector(nodeId))
  const isLoaded = useSelector(
    imageDataIsLoadedSelector(path ?? ''), // 応急処置
  )
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getImageData({ path, maxIndex }))
    }
  }, [dispatch, isLoaded, path, maxIndex])
  if (isLoaded && path != null) {
    return <ImagePlotImple path={path} />
  } else if (!isLoaded && path != null) {
    return <LinearProgress />
  } else {
    return null
  }
})

const ImagePlotImple = React.memo<{ path: string }>(({ path }) => {
  const activeIndex = useSelector(
    (state: RootState) => imageDataActiveIndexselector(path)(state) ?? 0,
  )
  const dispatch = useDispatch()
  const handleNext = () => dispatch(incrementImageActiveIndex({ path }))
  const handleBack = () => dispatch(decrementImageActiveIndex({ path }))
  const maxIndex = useSelector(imageDataMaxIndexSelector(path))
  const theme = useTheme()
  return (
    <>
      <MobileStepper
        steps={(maxIndex ?? 0) + 1}
        position="static"
        variant="text"
        activeStep={activeIndex}
        nextButton={
          <Button
            size="small"
            onClick={handleNext}
            disabled={activeIndex === (maxIndex ?? 0)}
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
      <ImagePlotChart path={path} />
    </>
  )
})

const ImagePlotChart = React.memo<{ path: string }>(({ path }) => {
  const imageData = useSelector(
    activeImageDataSelector(path),
    imageDataEqualtyFn,
  )
  const data = React.useMemo(
    () => [
      {
        z: imageData,
        type: 'heatmap',
        colorscale: [
          [0, '#000000'],
          [1, '#ffffff'],
        ],
        hoverongaps: false,
        showscale: false,
      },
    ],
    [imageData],
  )
  const layout = {
    title: 'dummy',
    margin: {
      t: 30, // top
      l: 90, // left
      b: 30, // bottom
    },
    autosize: true,
    height: 350,
    yaxis: {
      autorange: 'reversed',
    },
  }
  const config = {
    displayModeBar: true,
  }
  return <PlotlyChart data={data} layout={layout} config={config} />
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

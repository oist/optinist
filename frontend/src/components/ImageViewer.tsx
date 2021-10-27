import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import MobileStepper from '@material-ui/core/MobileStepper'
import Button from '@material-ui/core/Button'
import KeyboardArrowLeft from '@material-ui/icons/KeyboardArrowLeft'
import KeyboardArrowRight from '@material-ui/icons/KeyboardArrowRight'
import {
  currentImageMaxIndexSelector,
  currentImagePageIndexSelector,
  currentImageFolderSelector,
  currentImageFileNameSelector,
  currentImageIndexSelector,
  currentImageIsFulfilledSelector,
  currentImageBrightnessSelector,
  currentImageContrastSelector,
} from 'redux/slice/ImageIndex/ImageIndexSelector'
import Popover from '@material-ui/core/Popover'
import MoreVertIcon from '@material-ui/icons/MoreVert'
import Grid from '@material-ui/core/Grid'
import { RootState } from 'redux/store'
import {
  decrementPageIndex,
  incrementPageIndex,
  setBrightness,
  setContrast,
} from 'redux/slice/ImageIndex/ImageIndex'
import { LinearProgress, Slider, Typography, useTheme } from '@material-ui/core'
import { NodeIdContext } from 'App'

// import logo from './logo.svg';

export const ImageViewer = React.memo(() => {
  const nodeId = React.useContext(NodeIdContext)
  const currentImageIsSeleted = useSelector(
    (state: RootState) => currentImageIndexSelector(nodeId)(state) != null,
  )
  if (currentImageIsSeleted) {
    return <Viewer nodeId={nodeId} />
  } else {
    return null
  }
})

const Viewer = React.memo<{ nodeId: string }>(({ nodeId }) => {
  const dispatch = useDispatch()
  const maxIndex = useSelector(currentImageMaxIndexSelector(nodeId))
  const pageIndex = useSelector(currentImagePageIndexSelector(nodeId))
  const brightness = useSelector(currentImageBrightnessSelector(nodeId))
  const contrast = useSelector(currentImageContrastSelector(nodeId))
  const folder = useSelector(currentImageFolderSelector(nodeId))
  const fileName = useSelector(currentImageFileNameSelector(nodeId))
  const handleNext = () => dispatch(incrementPageIndex({ id: nodeId }))
  const handleBack = () => dispatch(decrementPageIndex({ id: nodeId }))
  const isLoaded = useSelector(currentImageIsFulfilledSelector(nodeId))
  const theme = useTheme()
  if (!isLoaded) {
    return <LinearProgress />
  }
  return (
    <div
      style={{
        height: '100%',
      }}
    >
      <MobileStepper
        steps={maxIndex ?? 0}
        position="static"
        variant="text"
        activeStep={pageIndex}
        nextButton={
          <Button
            size="small"
            onClick={handleNext}
            disabled={pageIndex === (maxIndex ?? 0) - 1}
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
          <Button size="small" onClick={handleBack} disabled={pageIndex === 0}>
            {theme.direction === 'rtl' ? (
              <KeyboardArrowRight />
            ) : (
              <KeyboardArrowLeft />
            )}
            <Typography>Back</Typography>
          </Button>
        }
      />
      <Grid
        container
        direction="row"
        justifyContent="space-between"
        alignItems="center"
      >
        <Grid item xs={2}>
          <OptionMenu />
        </Grid>
        <Grid item xs={8}>
          <Typography style={{ textAlign: 'center', alignItems: 'center' }}>
            {fileName}({nodeId})
          </Typography>
        </Grid>
        <Grid item xs={2} />
      </Grid>
      <div
        style={{
          textAlign: 'center',
          height: '80%',
        }}
      >
        <img
          style={{
            textAlign: 'center',
            height: '100%',
            maxWidth: '100%',
            filter: `contrast(${contrast}%) brightness(${brightness}%)`,
          }}
          alt=""
          src={`http://localhost:8000/api/${folder}/${pageIndex}.png`}
        />
      </div>
    </div>
  )
})

const OptionMenu = React.memo(() => {
  const nodeId = React.useContext(NodeIdContext)
  const dispatch = useDispatch()
  const menuAnchorEl = React.useRef<HTMLButtonElement>(null)
  const [open, setOpen] = React.useState(false)
  const handleClose = () => {
    setOpen(false)
  }
  const handleClick = () => {
    setOpen(true)
  }
  const brightness = useSelector(currentImageBrightnessSelector(nodeId))
  const contrast = useSelector(currentImageContrastSelector(nodeId))
  return (
    <>
      <Button
        size="small"
        variant="outlined"
        ref={menuAnchorEl}
        onClick={handleClick}
      >
        filter
      </Button>
      <Popover
        open={open}
        anchorEl={menuAnchorEl.current}
        onClose={handleClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'left',
        }}
        transformOrigin={{
          vertical: 'top',
          horizontal: 'center',
        }}
      >
        <div style={{ margin: 24, width: 160 }}>
          <div>
            <Typography>Brightness</Typography>
            <Slider
              value={brightness}
              onChange={(event, newValue) => {
                if (typeof newValue === 'number') {
                  dispatch(setBrightness({ id: nodeId, brightness: newValue }))
                }
              }}
              valueLabelDisplay="auto"
              marks={marks}
              max={300}
              min={0}
            />
          </div>
          <div style={{ marginTop: 16 }}>
            <Typography>Contrast</Typography>
            <Slider
              value={contrast}
              onChange={(event, newValue) => {
                if (typeof newValue === 'number') {
                  dispatch(setContrast({ id: nodeId, contrast: newValue }))
                }
              }}
              valueLabelDisplay="auto"
              marks={marks}
              max={300}
              min={0}
            />
          </div>
        </div>
      </Popover>
    </>
  )
})

const marks = [
  {
    value: 0,
    label: '0%',
  },
  {
    value: 100,
    label: '100%',
  },
  {
    value: 200,
    label: '200%',
  },
  {
    value: 300,
    label: '300%',
  },
]

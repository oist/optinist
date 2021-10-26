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
} from 'redux/slice/ImageIndex/ImageIndexSelector'

import { RootState } from 'redux/store'
import {
  decrementPageIndex,
  incrementPageIndex,
} from 'redux/slice/ImageIndex/ImageIndex'
import { LinearProgress, Typography, useTheme } from '@material-ui/core'
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
  const disptach = useDispatch()
  const maxIndex = useSelector(currentImageMaxIndexSelector(nodeId))
  const pageIndex = useSelector(currentImagePageIndexSelector(nodeId))
  const folder = useSelector(currentImageFolderSelector(nodeId))
  const fileName = useSelector(currentImageFileNameSelector(nodeId))
  const handleNext = () => disptach(incrementPageIndex({ id: nodeId }))
  const handleBack = () => disptach(decrementPageIndex({ id: nodeId }))
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
      <Typography style={{ textAlign: 'center' }}>
        {fileName}({nodeId})
      </Typography>
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
          }}
          alt=""
          src={`http://localhost:8000/api/${folder}/${pageIndex}.png`}
        />
      </div>
    </div>
  )
})

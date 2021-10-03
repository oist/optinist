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
  currentImageIdSelector,
} from 'redux/slice/ImageIndex/ImageIndexSelector'

import { RootState } from 'redux/store'
import {
  decrementPageIndex,
  incrementPageIndex,
} from 'redux/slice/ImageIndex/ImageIndex'
import { Typography, useTheme } from '@material-ui/core'
// import logo from './logo.svg';

export const ImageViewer = React.memo(() => {
  const currentImage = useSelector((state: RootState) =>
    currentImageIndexSelector(state),
  )
  if (currentImage !== undefined) {
    return <Viewer />
  } else {
    return null
  }
})

const Viewer = React.memo(() => {
  const disptach = useDispatch()
  const id = useSelector(currentImageIdSelector)
  const maxIndex = useSelector(currentImageMaxIndexSelector)
  const pageIndex = useSelector(currentImagePageIndexSelector)
  const folder = useSelector(currentImageFolderSelector)
  const fileName = useSelector(currentImageFileNameSelector)
  const handleNext = () => disptach(incrementPageIndex({ id }))
  const handleBack = () => disptach(decrementPageIndex({ id }))
  const theme = useTheme()
  return (
    <div style={{ padding: 8 }}>
      <Typography>
        {fileName}({id})
      </Typography>
      <img
        height={'100%'}
        src={`http://localhost:8000/files/${folder}/${pageIndex}.png`}
      />
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
    </div>
  )
})

import { nanoid } from '@reduxjs/toolkit'
import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { setFileUploaderStateById } from './FileUploaderSlice'
import { uploadFile } from './FileUploaderActions'
import {
  selectFileUploadIsPending,
  selectUploadFilePath,
  selectFileUploadIsFulfilled,
  selectFileUploadProgress,
  selectFileUploadIsUninitialized,
  selectFileUploadError,
} from './FileUploaderSelectors'
import { FILE_TYPE } from '../InputNode/InputNodeType'

export function useFileUploader(fileType?: FILE_TYPE) {
  const dispatch = useDispatch()
  const id = React.useRef(nanoid())
  React.useEffect(() => {
    dispatch(setFileUploaderStateById(id.current))
  }, [dispatch])
  const onUploadFile = React.useCallback(
    (formData: FormData, fileName: string) => {
      dispatch(
        uploadFile({ requestId: id.current, fileName, formData, fileType }),
      )
    },
    [dispatch, fileType],
  )
  const uninitialized = useSelector(selectFileUploadIsUninitialized(id.current))
  const filePath = useSelector(selectUploadFilePath(id.current))
  const pending = useSelector(selectFileUploadIsPending(id.current))
  const fulfilled = useSelector(selectFileUploadIsFulfilled(id.current))
  const progress = useSelector(selectFileUploadProgress(id.current))
  const error = useSelector(selectFileUploadError(id.current))
  return {
    filePath,
    uninitialized,
    pending,
    fulfilled,
    progress,
    error,
    onUploadFile,
  }
}

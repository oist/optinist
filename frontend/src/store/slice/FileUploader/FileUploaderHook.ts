import { nanoid } from '@reduxjs/toolkit'
import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
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
import { selectCurrentWorkspaceId } from '../Workspace/WorkspaceSelector'

type UseFileUploaderProps = {
  fileType?: FILE_TYPE
  nodeId?: string
}

export function useFileUploader({ fileType, nodeId }: UseFileUploaderProps) {
  const dispatch = useDispatch()
  const id = React.useRef(nanoid())
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const onUploadFile = React.useCallback(
    (formData: FormData, fileName: string) => {
      if (workspaceId) {
        dispatch(
          uploadFile({
            workspaceId,
            requestId: id.current,
            nodeId,
            fileName,
            formData,
            fileType,
          }),
        )
      } else {
        throw new Error('workspaceId is undefined')
      }
    },
    [dispatch, workspaceId, fileType, nodeId],
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

import { useCallback, useRef } from "react"
import { useDispatch, useSelector } from "react-redux"

import { nanoid } from "@reduxjs/toolkit"

import { uploadFile } from "store/slice/FileUploader/FileUploaderActions"
import {
  selectFileUploadIsPending,
  selectUploadFilePath,
  selectFileUploadIsFulfilled,
  selectFileUploadProgress,
  selectFileUploadIsUninitialized,
  selectFileUploadError,
} from "store/slice/FileUploader/FileUploaderSelectors"
import { FILE_TYPE } from "store/slice/InputNode/InputNodeType"
import { selectCurrentWorkspaceId } from "store/slice/Workspace/WorkspaceSelector"
import { AppDispatch } from "store/store"

type UseFileUploaderProps = {
  fileType?: FILE_TYPE
  nodeId?: string
}

export function useFileUploader({ fileType, nodeId }: UseFileUploaderProps) {
  const dispatch = useDispatch<AppDispatch>()
  const id = useRef(nanoid())
  const workspaceId = useSelector(selectCurrentWorkspaceId)
  const onUploadFile = useCallback(
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
        throw new Error("workspaceId is undefined")
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

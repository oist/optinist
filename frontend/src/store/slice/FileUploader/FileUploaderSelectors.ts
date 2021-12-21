import { RootState } from 'store/store'
import { inistialUploaderState } from './FileUploaderInitlalState'

export const selectFileUploader = (id: string) => (state: RootState) => {
  if (Object.keys(state.fileUploader).includes(id)) {
    return state.fileUploader[id]
  } else {
    return inistialUploaderState
  }
}

export const selectUploadFilePath = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).path

export const selectUploadFileName = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).fileName

export const selectFileUploadIsUninitialized =
  (id: string) => (state: RootState) =>
    selectFileUploader(id)(state).isUninitialized

export const selectFileUploadIsPending = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).pending

export const selectFileUploadIsFulfilled = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).fulfilled

export const selectFileUploadProgress = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).uploadingProgress

export const selectFileUploadError = (id: string) => (state: RootState) =>
  selectFileUploader(id)(state).error

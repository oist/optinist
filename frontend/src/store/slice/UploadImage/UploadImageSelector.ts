import { RootState } from '../../store'

export const uploadImageSelector = (state: RootState) => state.uploadImage

export const uploadImageByIdSelector = (id: string) => (state: RootState) => {
  const imageIndex = uploadImageSelector(state)
  if (Object.keys(imageIndex).includes(id)) {
    return imageIndex[id]
  } else {
    return undefined
  }
}

export const imageJsonPathByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.jsonPath

export const imageFileNameByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.fileName

export const imageIsUploadedByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.isFulfilled

export const imageIsUploadingByIdSelector =
  (id: string) => (state: RootState) =>
    uploadImageByIdSelector(id)(state)?.isUploading

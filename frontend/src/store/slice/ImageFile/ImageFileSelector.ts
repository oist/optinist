import { RootState } from '../../store'

export const imageFileSelector = (state: RootState) => state.imageFile

export const uploadImageByIdSelector = (id: string) => (state: RootState) => {
  const imageIndex = imageFileSelector(state)
  if (Object.keys(imageIndex).includes(id)) {
    return imageIndex[id]
  } else {
    return undefined
  }
}

export const imageTiffPathByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.path

export const imageFileNameByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.fileName

export const imageMaxIndexByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.maxIndex

export const imageIsUploadedByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.isFulfilled

export const imageIsUploadingByIdSelector =
  (id: string) => (state: RootState) =>
    uploadImageByIdSelector(id)(state)?.isUploading

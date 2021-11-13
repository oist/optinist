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

export const imagePathByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.path

export const imageIsUploadedByIdSelector = (id: string) => (state: RootState) =>
  uploadImageByIdSelector(id)(state)?.isFulfilled

import { RootState } from '../../store'

export const imageFileSelector = (state: RootState) => state.fileData.image
export const csvFileSelector = (state: RootState) => state.fileData.csv

export const imageFileByIdSelector = (id: string) => (state: RootState) => {
  const imageFile = imageFileSelector(state)
  if (Object.keys(imageFile).includes(id)) {
    return imageFile[id]
  } else {
    return undefined
  }
}

export const csvFileByIdSelector = (id: string) => (state: RootState) => {
  const csvFile = csvFileSelector(state)
  if (Object.keys(csvFile).includes(id)) {
    return csvFile[id]
  } else {
    return undefined
  }
}

export const imageTiffPathByIdSelector = (id: string) => (state: RootState) =>
  imageFileByIdSelector(id)(state)?.path

export const imageFileNameByIdSelector = (id: string) => (state: RootState) =>
  imageFileByIdSelector(id)(state)?.fileName

export const imageMaxIndexByIdSelector = (id: string) => (state: RootState) =>
  imageFileByIdSelector(id)(state)?.maxIndex

export const imageIsUploadedByIdSelector = (id: string) => (state: RootState) =>
  imageFileByIdSelector(id)(state)?.isFulfilled

export const imageIsUploadingByIdSelector =
  (id: string) => (state: RootState) =>
    imageFileByIdSelector(id)(state)?.isUploading

export const csvIsUploadedByIdSelector = (id: string) => (state: RootState) =>
  csvFileByIdSelector(id)(state)?.isFulfilled

export const csvIsUploadingByIdSelector = (id: string) => (state: RootState) =>
  csvFileByIdSelector(id)(state)?.isUploading

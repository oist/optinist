import { RootState } from '../../store'

export const imageIndexSelector = (state: RootState) => state.imageIndex

// export const currentImageIdSelector = (state: RootState) =>
//   imageIndexSelector(state).currentImageId

export const imageIndexByIdSelector = (id: string) => (state: RootState) => {
  const imageIndex = imageIndexSelector(state)
  if (Object.keys(imageIndex.index).includes(id)) {
    return imageIndex.index[id]
  } else {
    return undefined
  }
}

export const currentImageIndexSelector = (id: string) => (state: RootState) => {
  try {
    return imageIndexByIdSelector(id)(state)
  } catch (e) {
    return undefined
  }
}

export const currentImageMaxIndexSelector =
  (id: string) => (state: RootState) =>
    currentImageIndexSelector(id)(state)?.maxIndex

export const currentImageFolderSelector = (id: string) => (state: RootState) =>
  currentImageIndexSelector(id)(state)?.folder

export const currentImageFileNameSelector =
  (id: string) => (state: RootState) =>
    currentImageIndexSelector(id)(state)?.fileName

export const currentImagePageIndexSelector =
  (id: string) => (state: RootState) =>
    currentImageIndexSelector(id)(state)?.pageIndex

export const currentImageIsFulfilledSelector =
  (id: string) => (state: RootState) =>
    currentImageIndexSelector(id)(state)?.isFulfilled

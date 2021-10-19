import { RootState } from '../../store'

export const imageIndexSelector = (state: RootState) => state.imageIndex

export const currentImageIdSelector = (state: RootState) =>
  imageIndexSelector(state).currentImageId

export const imageIndexByIdSelector = (id: string) => (state: RootState) => {
  const imageIndex = imageIndexSelector(state)
  if (Object.keys(imageIndex.index).includes(id)) {
    return imageIndex.index[id]
  } else {
    return undefined
  }
}

export const currentImageIndexSelector = (state: RootState) => {
  try {
    console.log(state)
    const currentNodeId = currentImageIdSelector(state)
    return imageIndexByIdSelector(currentNodeId)(state)
  } catch (e) {
    return undefined
  }
}

export const currentImageMaxIndexSelector = (state: RootState) =>
  currentImageIndexSelector(state)?.maxIndex

export const currentImageFolderSelector = (state: RootState) =>
  currentImageIndexSelector(state)?.folder

export const currentImageFileNameSelector = (state: RootState) =>
  currentImageIndexSelector(state)?.fileName

export const currentImagePageIndexSelector = (state: RootState) =>
  currentImageIndexSelector(state)?.pageIndex

export const currentImageIsFulfilledSelector = (state: RootState) =>
  currentImageIndexSelector(state)?.isFulfilled

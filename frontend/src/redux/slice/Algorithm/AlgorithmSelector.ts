import { RootState } from '../../store'

export const algorithmSelector = (state: RootState) => state.algorithm

export const currentAlgoIdSelector = (state: RootState) =>
  state.algorithm.currentAlgoId

export const currentAlgoNameSelector = (state: RootState) => {
  const id = currentAlgoIdSelector(state)
  if (Object.keys(state.algorithm.algoMap).includes(id)) {
    return state.algorithm.algoMap[id].name
  } else {
    return undefined
  }
}

export const algoParamByIdSelector = (id: string) => (state: RootState) => {
  const algoMap = algorithmSelector(state).algoMap
  if (Object.keys(algoMap).includes(id)) {
    return algoMap[id].param
  } else {
    return undefined
  }
}

export const paramValueSelector =
  (id: string, paramName: string) => (state: RootState) => {
    const param = algorithmSelector(state).algoMap[id].param
    if (param != null) {
      return param[paramName]
    } else {
      return null
    }
  }

export const outputPathListSelector = (id: string) => (state: RootState) => {
  if (Object.keys(state.algorithm.algoMap).includes(id)) {
    const outputPaths = state.algorithm.algoMap[id].output
    if (outputPaths != null) {
      return Object.entries(outputPaths)
    } else {
      return []
    }
  } else {
    return []
  }
}

export const selectedOutputPathSelector =
  (id: string) => (state: RootState) => {
    if (Object.keys(state.algorithm.algoMap).includes(id)) {
      return state.algorithm.algoMap[id].selectedPath
    } else {
      return null
    }
  }

export const imageDirMaxIndexByIdSelector =
  (id: string) => (state: RootState) => {
    if (Object.keys(state.algorithm.algoMap).includes(id)) {
      return state.algorithm.algoMap[id].output?.images?.maxIndex ?? null
    } else {
      return null
    }
  }

export const currentOutputDataSelector = (state: RootState) => {
  const id = currentAlgoIdSelector(state)
  if (Object.keys(state.algorithm.plotDataMap).includes(id)) {
    return state.algorithm.plotDataMap[id]
  } else {
    return undefined
  }
}

export const outputDataIsLoadedByIdSelector =
  (id: string) => (state: RootState) => {
    return Object.keys(state.algorithm.plotDataMap).includes(id)
  }

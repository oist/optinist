import { RootState } from '../../store'
import { isImageOutput } from './AlgorithmUtils'

export const algorithmSelector = (state: RootState) => state.algorithm

export const currentAlgoIdSelector = (state: RootState) =>
  state.algorithm.currentAlgoId

export const algoNameByIdSelector = (id: string) => (state: RootState) => {
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

export const outputKeyListSelector = (id: string) => (state: RootState) => {
  if (Object.keys(state.algorithm.algoMap).includes(id)) {
    const outputPaths = state.algorithm.algoMap[id].output
    if (outputPaths != null) {
      return Object.keys(outputPaths)
    } else {
      return []
    }
  } else {
    return []
  }
}

export const selectedOutputKeySelector = (id: string) => (state: RootState) => {
  if (Object.keys(state.algorithm.algoMap).includes(id)) {
    return state.algorithm.algoMap[id].selectedOutputKey
  } else {
    return undefined
  }
}

const selectedOutputPathSelector = (id: string) => (state: RootState) => {
  const algoMap = state.algorithm.algoMap
  if (algoMap[id] != null) {
    const algo = algoMap[id]
    const output = algo.output
    const selectedKey = selectedOutputKeySelector(id)(state)
    if (output != null && selectedKey != null) {
      return output[selectedKey]
    }
  }
  return undefined
}

export const selectedOutputPathTypeSelector =
  (id: string) => (state: RootState) => {
    return selectedOutputPathSelector(id)(state)?.type
  }

export const selectedOutputPathValueSelector =
  (id: string) => (state: RootState) => {
    return selectedOutputPathSelector(id)(state)?.path.value
  }

export const imagePathMaxIndexByIdSelector =
  (id: string, outputKey: string) => (state: RootState) => {
    if (Object.keys(state.algorithm.algoMap).includes(id)) {
      const outputPaths = state.algorithm.algoMap[id].output
      if (outputPaths != null && Object.keys(outputPaths).includes(outputKey)) {
        const path = outputPaths[outputKey]
        if (isImageOutput(path)) {
          return path.path.maxIndex
        } else {
          return null
        }
      } else {
        return null
      }
    } else {
      return null
    }
  }

export const outputPathByIdSelector =
  (nodeId: string, outputKey: string) => (state: RootState) => {
    const output = state.algorithm.algoMap[nodeId]
    if (output != null && output.output != null) {
      const outputPath = output.output[outputKey]
      if (outputPath != null) {
        return outputPath
      }
    }
    return null
  }

export const outputPathTypeByIdSelector =
  (nodeId: string, outputKey: string) => (state: RootState) => {
    const outputPath = outputPathByIdSelector(nodeId, outputKey)(state)
    if (outputPath != null) {
      return outputPath.type
    } else {
      return null
    }
  }

export const outputPathValueByIdSelector =
  (nodeId: string, outputKey: string) => (state: RootState) => {
    const outputPath = outputPathByIdSelector(nodeId, outputKey)(state)
    if (outputPath != null) {
      return outputPath.path.value
    }
    return null
  }

import { RootState } from 'store/store'
import { selectAlgorithmName } from '../AlgorithmNode/AlgorithmNodeSelectors'
import { getAlgoChild } from './AlgorithmListUtils'

export const selectAlgorithmList = (state: RootState) => state.algorithmList

export const selectAlgorithmListIsLated = (state: RootState) =>
  selectAlgorithmList(state).isLatest

export const selectAlgorithmListTree = (state: RootState) =>
  selectAlgorithmList(state).tree

export const selectAlgoArgs = (nodeId: string) => (state: RootState) => {
  const algoName = selectAlgorithmName(nodeId)(state)
  if (algoName != null) {
    const algoListChild = getAlgoChild(selectAlgorithmListTree(state), algoName)
    return algoListChild?.args
  } else {
    return undefined
  }
}

export const selectAlgoReturns = (nodeId: string) => (state: RootState) => {
  const algoName = selectAlgorithmName(nodeId)(state)
  if (algoName != null) {
    const algoListChild = getAlgoChild(selectAlgorithmListTree(state), algoName)
    return algoListChild?.returns
  } else {
    return undefined
  }
}

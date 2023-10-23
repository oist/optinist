import { getAlgoChild } from "store/slice/AlgorithmList/AlgorithmListUtils"
import { selectAlgorithmName } from "store/slice/AlgorithmNode/AlgorithmNodeSelectors"
import { RootState } from "store/store"

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

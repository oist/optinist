import { RootState } from 'store/store'

import { getChildParam } from './AlgorithmNodeUtils'

export const selectAlgorithmNode = (state: RootState) => state.algorithmNode

export const selectAlgorithmNodeDefined =
  (nodeId: string) => (state: RootState) =>
    Object.keys(state.algorithmNode).includes(nodeId)

export const selectAlgorithmFunctionPath =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNode(state)[nodeId].functionPath

export const selectAlgorithmName = (nodeId: string) => (state: RootState) =>
  selectAlgorithmNode(state)[nodeId].name

export const selectAlgorithmParams = (nodeId: string) => (state: RootState) =>
  selectAlgorithmNode(state)[nodeId].params

export const selectAlgorithmParamsExit =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNode(state)[nodeId].params !== null

export const selectAlgorithmParamsKeyList =
  (nodeId: string) => (state: RootState) =>
    Object.keys(selectAlgorithmNode(state)[nodeId]?.params ?? {})

export const selectAlgorithmParamsValue =
  (nodeId: string, path: string) => (state: RootState) => {
    const params = selectAlgorithmParams(nodeId)(state)
    if (params != null) {
      const target = getChildParam(path, params)
      return target?.value
    } else {
      throw new Error('AlgorithmParam is null')
    }
  }

export const selectAlgorithmParam =
  (nodeId: string, paramKey: string) => (state: RootState) => {
    const params = selectAlgorithmParams(nodeId)(state)
    if (params != null) {
      return params[paramKey]
    } else {
      throw new Error('AlgorithmParam is null')
    }
  }

export const selectAlgorithmSelectedOutputKey =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNode(state)[nodeId].selectedOutputKey

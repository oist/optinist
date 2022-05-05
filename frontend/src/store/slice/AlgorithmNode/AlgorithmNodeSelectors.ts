import { RootState } from 'store/store'

import { getChildParam } from 'store/utils/param/ParamUtils'

export const selectAlgorithmNode = (state: RootState) => state.algorithmNode

export const selectAlgorithmNodeById = (nodeId: string) => (state: RootState) =>
  state.algorithmNode[nodeId]

export const selectAlgorithmNodeDefined =
  (nodeId: string) => (state: RootState) =>
    Object.keys(state.algorithmNode).includes(nodeId)

export const selectAlgorithmFunctionPath =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNodeById(nodeId)(state).functionPath

export const selectAlgorithmName = (nodeId: string) => (state: RootState) =>
  selectAlgorithmNodeById(nodeId)(state).name

export const selectAlgorithmParams = (nodeId: string) => (state: RootState) =>
  selectAlgorithmNodeById(nodeId)(state).params

export const selectAlgorithmIsUpdated =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNodeById(nodeId)(state).isUpdated

export const selectAlgorithmParamsExit =
  (nodeId: string) => (state: RootState) =>
    selectAlgorithmNodeById(nodeId)(state).params !== null

export const selectAlgorithmParamsKeyList =
  (nodeId: string) => (state: RootState) =>
    Object.keys(selectAlgorithmNodeById(nodeId)(state)?.params ?? {})

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

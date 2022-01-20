import { RootState } from 'store/store'
import { selectAlgorithmName } from '../AlgorithmNode/AlgorithmNodeSelectors'

export const selectRunPipelineResult = (state: RootState) =>
  state.runPipelineResult.result

export const selectOutputPaths = (state: RootState) => {
  const result = selectRunPipelineResult(state)
  if (
    result.outputPaths != null &&
    Object.keys(result.outputPaths).length > 0
  ) {
    return result.outputPaths
  } else {
    return null
  }
}

export const selectResultError = (nodeId: string) => (state: RootState) => {
  const algoName = selectAlgorithmName(nodeId)(state)
  const result = selectRunPipelineResult(state)
  if (!!result.name && result.name === algoName) {
    return result?.message
  } else {
    return null
  }
}

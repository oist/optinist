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
  (id: string, paramName: string) => (state: RootState) =>
    algorithmSelector(state).algoMap[id].param[paramName]

// export const currentOutputDataSelector = (state: RootState) => {
//   const id = currentAlgoIdSelector(state)
//   if (Object.keys(state.algorithm.algoMap).includes(id)) {
//     return state.algorithm.algoMap[id].output?.data
//   } else {
//     return undefined
//   }
// }

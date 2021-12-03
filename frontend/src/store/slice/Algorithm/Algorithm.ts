import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { isAlgoNodeData } from 'utils/ElementUtils'
import { addFlowElement } from '../Element/Element'
import {
  getAlgoList,
  getAlgoParams,
  reflectRunPipelineResult,
} from './AlgorithmAction'
import {
  ALGORITHM_SLICE_NAME,
  Algorithm,
  OUTPUT_TYPE_SET,
} from './AlgorithmType'
import { convertToAlgoListType } from './AlgorithmUtils'

const initialState: Algorithm = {
  algoNodeMap: {},
  algoList: {},
}

export const algorithmSlice = createSlice({
  name: ALGORITHM_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        id: string
        paramKey: string
        newValue: unknown
      }>,
    ) => {
      const { id, paramKey, newValue } = action.payload
      const param = state.algoNodeMap[id].param
      if (param !== undefined) {
        param[paramKey] = newValue
      }
    },
    setSelectedOutputKey: (
      state,
      action: PayloadAction<{
        id: string
        outputKey: string
      }>,
    ) => {
      state.algoNodeMap[action.payload.id].selectedOutputKey =
        action.payload.outputKey
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(addFlowElement, (state, action) => {
        if (isAlgoNodeData(action.payload)) {
          const algoName = action.payload.data?.label
          const path = action.payload.data?.path
          if (algoName != null && path != null) {
            state.algoNodeMap[action.payload.id] = {
              name: algoName,
              path,
            }
          }
        }
      })
      .addCase(getAlgoList.fulfilled, (state, action) => {
        state.algoList = convertToAlgoListType(action.payload)
      })
      .addCase(getAlgoParams.fulfilled, (state, action) => {
        const { id, algoName } = action.meta.arg
        const algo = state.algoNodeMap[id]
        state.algoNodeMap[id] = {
          ...algo,
          name: algoName,
          param: action.payload,
        }
      })
      .addCase(reflectRunPipelineResult, (state, action) => {
        Object.entries(action.payload.dto).forEach(([path, outputPaths]) => {
          Object.entries(state.algoNodeMap).forEach(([id, algo]) => {
            // todo とりあえず名前一致だが、後でサーバーサイドとフロントで両方idにする
            if (algo.path === path && state.algoNodeMap[id]) {
              const outputState = {
                ...state.algoNodeMap[id].output,
              }
              Object.entries(outputPaths).forEach(([key, pathInfo]) => {
                if (pathInfo.type === 'images') {
                  outputState[key] = {
                    type: OUTPUT_TYPE_SET.IMAGE,
                    path: {
                      value: pathInfo.path,
                    },
                  }
                } else if (pathInfo.type === 'timeseries') {
                  outputState[key] = {
                    type: OUTPUT_TYPE_SET.TIME_SERIES,
                    path: {
                      value: pathInfo.path,
                    },
                  }
                } else if (pathInfo.type === 'heatmap') {
                  outputState[key] = {
                    type: OUTPUT_TYPE_SET.HEAT_MAP,
                    path: {
                      value: pathInfo.path,
                    },
                  }
                }
                state.algoNodeMap[id].selectedOutputKey = key
              })
              state.algoNodeMap[id].output = outputState
            }
          })
        })
        Object.entries(state.algoNodeMap).forEach(([id, algo]) => {
          if (algo.name === action.payload.error?.name) {
            state.algoNodeMap[id].error = action.payload.error?.message
          } else {
            state.algoNodeMap[id].error = undefined
          }
        })
      })
  },
})

export const { updateParam, setSelectedOutputKey } = algorithmSlice.actions

export default algorithmSlice.reducer

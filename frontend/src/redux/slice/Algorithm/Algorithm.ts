import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData } from 'utils/ElementUtils'
import { addFlowElement } from '../Element/Element'
import { clickNode, runPipeline } from '../Element/ElementAction'
import { getAlgoList, getAlgoParams } from './AlgorithmAction'
import {
  ALGORITHM_SLICE_NAME,
  Algorithm,
  OUTPUT_TYPE_SET,
} from './AlgorithmType'

const initialState: Algorithm = {
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoNodeMap: {},
  algoList: {},
}

export const algorithmSlice = createSlice({
  name: ALGORITHM_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{ paramKey: string; newValue: unknown }>,
    ) => {
      const { paramKey, newValue } = action.payload
      const param = state.algoNodeMap[state.currentAlgoId].param
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
      .addCase(clickNode, (state, action) => {
        if (action.payload.type === NODE_DATA_TYPE_SET.ALGO) {
          state.currentAlgoId = action.payload.id
        }
      })
      .addCase(addFlowElement, (state, action) => {
        if (isAlgoNodeData(action.payload)) {
          state.algoNodeMap[action.payload.id] = {
            name: action.payload.data?.label ?? '',
          }
        }
      })
      .addCase(getAlgoList.fulfilled, (state, action) => {
        state.algoList = action.payload
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
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          Object.entries(action.payload.outputPaths).forEach(
            ([algoName, outputPaths]) => {
              // console.log(current(state.algoMap))
              Object.entries(state.algoNodeMap).forEach(([id, algo]) => {
                // todo とりあえず名前一致だが、後でサーバーサイドとフロントで両方idにする
                if (algo.name === algoName && state.algoNodeMap[id]) {
                  const outputState = {
                    ...state.algoNodeMap[id].output,
                  }
                  Object.entries(outputPaths).forEach(([key, pathInfo]) => {
                    if (pathInfo.type === 'images') {
                      outputState[key] = {
                        type: OUTPUT_TYPE_SET.IMAGE,
                        path: {
                          value: pathInfo.path,
                          maxIndex: pathInfo.max_index ?? 0,
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
            },
          )
        }
      })
  },
})

export const { updateParam, setSelectedOutputKey } = algorithmSlice.actions

export default algorithmSlice.reducer

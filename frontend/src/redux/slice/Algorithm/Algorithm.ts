import { createSlice, current, PayloadAction } from '@reduxjs/toolkit'

import { INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData } from 'utils/ElementUtils'
import { addFlowElement } from '../Element/Element'
import { clickNode, runPipeline } from '../Element/ElementAction'
import { getAlgoOutputData, getAlgoParams } from './AlgorithmAction'
import { convertToOutputData } from './AlgorithmUtils'
import { ALGORITHM_SLICE_NAME, Algorithm } from './AlgorithmType'

const initialState: Algorithm = {
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoMap: {},
  plotDataMap: {},
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
      const param = state.algoMap[state.currentAlgoId].param
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
      state.algoMap[action.payload.id].selectedOutputKey =
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
          state.algoMap[action.payload.id] = {
            name: action.payload.data?.label ?? '',
          }
        }
      })
      .addCase(getAlgoParams.fulfilled, (state, action) => {
        const { id, algoName } = action.meta.arg
        const algo = state.algoMap[id]
        state.algoMap[id] = {
          ...algo,
          name: algoName,
          param: action.payload,
        }
      })
      .addCase(getAlgoOutputData.fulfilled, (state, action) => {
        state.plotDataMap[action.meta.arg.id] = convertToOutputData(
          action.payload,
        )
      })
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          Object.entries(action.payload.outputPaths).forEach(
            ([algoName, outputPaths]) => {
              // console.log(current(state.algoMap))
              Object.entries(state.algoMap).forEach(([id, algo]) => {
                // todo とりあえず名前一致だが、後でサーバーサイドとフロントで両方idにする
                if (algo.name === algoName && state.algoMap[id]) {
                  const outputState = {
                    ...state.algoMap[id].output,
                  }
                  Object.entries(outputPaths).forEach(([key, pathInfo]) => {
                    if (pathInfo.type === 'images') {
                      outputState[key] = {
                        type: 'image',
                        path: {
                          value: pathInfo.path,
                          maxIndex: pathInfo.max_index ?? 0,
                        },
                      }
                    } else if (pathInfo.type === 'timeseries') {
                      outputState[key] = {
                        type: 'plotData',
                        path: {
                          value: pathInfo.path,
                        },
                      }
                    }
                    state.algoMap[id].selectedOutputKey = key
                  })
                  state.algoMap[id].output = outputState
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

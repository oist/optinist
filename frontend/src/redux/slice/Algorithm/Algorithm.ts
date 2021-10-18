import { createSlice, PayloadAction } from '@reduxjs/toolkit'

import { INITIAL_ALGO_ELEMENT_ID } from 'const/flowchart'
import { NODE_DATA_TYPE_SET } from 'const/NodeData'
import { isAlgoNodeData } from 'utils/ElementUtils'
import { addFlowElement } from '../Element/Element'
import { clickNode, runPipeline } from '../Element/ElementAction'
import { getAlgoOutputData, getAlgoParams } from './AlgorithmAction'

import { ALGORITHM_SLICE_NAME, Algorithm } from './AlgorithmType'

const initialState: Algorithm = {
  currentAlgoId: INITIAL_ALGO_ELEMENT_ID,
  algoMap: {},
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
            // param: {},
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
      // .addCase(getAlgoOutputData.fulfilled, (state, action) => {
      //   console.log(action.payload)
      //   state.algoMap[action.meta.arg.id].output = { data: action.payload }
      // })
      .addCase(runPipeline.fulfilled, (state, action) => {
        if (action.payload.message === 'success') {
          Object.entries(action.payload.outputPaths).forEach(([name, dir]) => {
            Object.entries(state.algoMap).forEach(([id, algo]) => {
              // todo とりあえず名前一致だが、後でサーバーサイドとフロントで両方idにする
              if (algo.name === name && state.algoMap[id]) {
                state.algoMap[id].output = {
                  imageDir: dir.image_dir,
                }
              }
            })
          })
        }
      })
  },
})

export const { updateParam } = algorithmSlice.actions

export default algorithmSlice.reducer

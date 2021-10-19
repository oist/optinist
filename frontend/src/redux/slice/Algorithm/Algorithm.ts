import { createSlice, PayloadAction } from '@reduxjs/toolkit'

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
    setSelectedOutputPath: (
      state,
      action: PayloadAction<{
        id: string
        path: string
      }>,
    ) => {
      // todo もっといい状態の持ち方と判定方法があるはず...
      if (action.payload.path.includes('images')) {
        state.algoMap[action.payload.id].selectedPath = {
          value: action.payload.path,
          isImage: true,
        }
      } else {
        state.algoMap[action.payload.id].selectedPath = {
          value: action.payload.path,
          isImage: false,
        }
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
            selectedPath: null,
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
          Object.entries(action.payload.outputPaths).forEach(([name, dir]) => {
            Object.entries(state.algoMap).forEach(([id, algo]) => {
              // todo とりあえず名前一致だが、後でサーバーサイドとフロントで両方idにする
              if (algo.name === name && state.algoMap[id]) {
                state.algoMap[id].output = {}
                if (dir.image_dir != null) {
                  state.algoMap[id].output = {
                    images: {
                      path: dir.image_dir?.path ?? '',
                      maxIndex: dir.image_dir?.max_index ?? 0,
                    },
                  }
                }
                if (dir.fluo_path != null) {
                  state.algoMap[id].output = {
                    ...state.algoMap[id].output,
                    fluo: dir.fluo_path,
                  }
                }
                state.algoMap[id].selectedPath = {
                  value: dir.image_dir?.path ?? null,
                  isImage: true,
                } // 本来は意味のあるkeyを使用する
              }
            })
          })
        }
      })
  },
})

export const { updateParam, setSelectedOutputPath } = algorithmSlice.actions

export default algorithmSlice.reducer

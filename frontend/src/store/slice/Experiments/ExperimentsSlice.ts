import { createSlice, isAnyOf } from "@reduxjs/toolkit"

import {
  getExperiments,
  deleteExperimentByUid,
  deleteExperimentByList,
} from "store/slice/Experiments/ExperimentsActions"
import { EXPERIMENTS_SLICE_NAME, Experiments } from "store/slice/Experiments/ExperimentsType"
import {
  convertToExperimentListType,
  convertToExperimentType,
} from "store/slice/Experiments/ExperimentsUtils"
import {
  pollRunResult,
  run,
  runByCurrentUid,
} from "store/slice/Pipeline/PipelineActions"
import { fetchWorkflow, reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"


export const initialState: Experiments = {
  status: "uninitialized",
}

export const experimentsSlice = createSlice({
  name: EXPERIMENTS_SLICE_NAME,
  initialState: initialState as Experiments,
  reducers: {
    clearExperiments: () => initialState,
  },
  extraReducers: (builder) => {
    builder
      .addCase(getExperiments.pending, () => {
        return {
          status: "pending",
        }
      })
      .addCase(getExperiments.fulfilled, (state, action) => {
        const experimentList = convertToExperimentListType(action.payload)
        return {
          status: "fulfilled",
          experimentList,
        }
      })
      .addCase(getExperiments.rejected, (state, action) => {
        return {
          status: "error",
          message: action.error.message,
        }
      })
      .addCase(deleteExperimentByUid.fulfilled, (state, action) => {
        if (action.payload && state.status === "fulfilled") {
          delete state.experimentList[action.meta.arg]
        }
      })
      .addCase(deleteExperimentByList.fulfilled, (state, action) => {
        if (action.payload && state.status === "fulfilled") {
          action.meta.arg.map((v) => delete state.experimentList[v])
        }
      })
      .addCase(pollRunResult.fulfilled, (state, action) => {
        if (state.status === "fulfilled") {
          const uid = action.meta.arg.uid
          const target = state.experimentList[uid]
          Object.entries(action.payload).forEach(([nodeId, value]) => {
            if (value.status === "success") {
              target.functions[nodeId].status = "success"
            } else if (value.status === "error") {
              target.functions[nodeId].status = "error"
            }
          })
        }
      })
      .addMatcher(
        isAnyOf(fetchWorkflow.fulfilled, reproduceWorkflow.fulfilled),
        (state, action) => {
          if (state.status === "fulfilled") {
            state.experimentList[action.payload.unique_id] =
              convertToExperimentType(action.payload)
          }
        },
      )
      .addMatcher(isAnyOf(run.fulfilled, runByCurrentUid.fulfilled), () => {
        return {
          status: "uninitialized",
        }
      })
  },
})

export const { clearExperiments } = experimentsSlice.actions
export default experimentsSlice.reducer

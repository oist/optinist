import { createSlice, isAnyOf } from "@reduxjs/toolkit"

import {
  getExperiments,
  deleteExperimentByUid,
  deleteExperimentByList,
  copyExperimentByList,
} from "store/slice/Experiments/ExperimentsActions"
import {
  EXPERIMENTS_SLICE_NAME,
  Experiments,
} from "store/slice/Experiments/ExperimentsType"
import {
  convertToExperimentListType,
  convertToExperimentType,
} from "store/slice/Experiments/ExperimentsUtils"
import {
  pollRunResult,
  run,
  runByCurrentUid,
} from "store/slice/Pipeline/PipelineActions"
import {
  fetchWorkflow,
  importSampleData,
  reproduceWorkflow,
} from "store/slice/Workflow/WorkflowActions"

export const initialState: Experiments = {
  status: "uninitialized",
  loading: true,
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
          loading: true,
        }
      })
      .addCase(getExperiments.fulfilled, (state, action) => {
        const experimentList = convertToExperimentListType(action.payload)
        return {
          status: "fulfilled",
          experimentList,
          loading: false,
        }
      })
      .addCase(getExperiments.rejected, (state, action) => {
        return {
          status: "error",
          message: action.error.message,
          loading: false,
        }
      })
      .addCase(deleteExperimentByUid.fulfilled, (state, action) => {
        state.loading = false
        if (action.payload && state.status === "fulfilled") {
          delete state.experimentList[action.meta.arg]
        }
      })
      .addCase(deleteExperimentByList.fulfilled, (state, action) => {
        state.loading = false
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
      .addCase(importSampleData.pending, (state) => {
        state.loading = true
      })
      .addCase(importSampleData.fulfilled, (state, _) => {
        state.loading = false
      })
      .addCase(importSampleData.rejected, (state, action) => {
        return {
          status: "error",
          message: action.error.message,
          loading: false,
        }
      })
      .addMatcher(
        isAnyOf(
          deleteExperimentByUid.pending,
          deleteExperimentByList.pending,
          copyExperimentByList.pending,
        ),
        (state) => {
          state.loading = true
        },
      )
      .addMatcher(
        isAnyOf(
          deleteExperimentByUid.rejected,
          deleteExperimentByList.rejected,
          copyExperimentByList.fulfilled,
          copyExperimentByList.rejected,
        ),
        (state) => {
          state.loading = false
        },
      )
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

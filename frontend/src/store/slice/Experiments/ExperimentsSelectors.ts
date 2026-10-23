import { EXPERIMENTS_STATUS } from "store/slice/Experiments/ExperimentsType"
import { RootState } from "store/store"

export const selectExperiments = (state: RootState) => state.experiments

export const selectExperimentsStatusIsUninitialized = (state: RootState) =>
  selectExperiments(state).status === "uninitialized"

export const selectExperimentsStatusIsPending = (state: RootState) =>
  selectExperiments(state).status === "pending"

export const selectExperimentsStatusIsFulfilled = (state: RootState) =>
  selectExperiments(state).status === "fulfilled"

export const selectExperimentsStatusIsError = (state: RootState) =>
  selectExperiments(state).status === "error"

export const selectExperimentsErrorMessage = (state: RootState) => {
  const experiments = selectExperiments(state)
  if (experiments.status === "error") {
    return experiments.message
  } else {
    throw new Error("experiments status is not error")
  }
}

export const selectLoading = (state: RootState) =>
  selectExperiments(state).loading

export const selectExperimentList = (state: RootState) => {
  const experiments = selectExperiments(state)
  if (experiments.status === "fulfilled") {
    return experiments.experimentList
  } else {
    throw new Error("experiments status is not fulfilled")
  }
}

export const selectExperimentUidList = (state: RootState) =>
  Object.keys(selectExperimentList(state))

export const selectExperiment = (uid: string) => (state: RootState) =>
  selectExperimentList(state)[uid]

export const selectExperimentStartedAt = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).startedAt

export const selectExperimentDataUsage = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).data_usage

export const selectExperimentFinishedAt = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).finishedAt

export const selectExperimentName = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).name

export const selectExperimentHasNWB = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).hasNWB

export const selectExperimentStatus =
  (uid: string) =>
  (state: RootState): EXPERIMENTS_STATUS => {
    const experiment = selectExperimentList(state)[uid]
    if (experiment.status) {
      return experiment.status
    }

    const functions = selectExperimentList(state)[uid].functions
    const statusList = Object.values(functions).map((f) => f.status)
    if (statusList.findIndex((status) => status === "error") >= 0) {
      return "error"
    } else if (statusList.findIndex((status) => status === "running") >= 0) {
      return "running"
    } else {
      return "success"
    }
  }

export const selectExperimentCheckList =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentFunction(uid, nodeId)(state).status

export const selectExperimentFunctionNodeIdList =
  (uid: string) => (state: RootState) =>
    Object.keys(selectExperimentList(state)[uid].functions)

export const selectExperimentFunction =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentList(state)[uid].functions[nodeId]

export const selectExperimentFunctionName =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentFunction(uid, nodeId)(state).name

export const selectExperimentFunctionStatus =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentFunction(uid, nodeId)(state).status

export const selectExperimentFunctionMessage =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentFunction(uid, nodeId)(state).message

export const selectExperimentFunctionHasNWB =
  (uid: string, nodeId: string) => (state: RootState) =>
    selectExperimentFunction(uid, nodeId)(state).hasNWB

export const selectFrameRate =
  (currentPipelineUid?: string) => (state: RootState) => {
    // Check if experiments have been acquired
    if (currentPipelineUid && selectExperimentsStatusIsFulfilled(state)) {
      return selectExperiment(currentPipelineUid)(state).frameRate || 50
    } else {
      return 50
    }
  }

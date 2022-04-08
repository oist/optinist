import { RootState } from 'store/store'
import { EXPERIMENTS_STATUS } from './ExperimentsType'

const selectExperiments = (state: RootState) => state.experiments

export const selectExperimentsSatusIsUninitialized = (state: RootState) =>
  selectExperiments(state).status === 'uninitialized'

export const selectExperimentsSatusIsPending = (state: RootState) =>
  selectExperiments(state).status === 'pending'

export const selectExperimentsSatusIsFulfilled = (state: RootState) =>
  selectExperiments(state).status === 'fulfilled'

export const selectExperimentsSatusIsError = (state: RootState) =>
  selectExperiments(state).status === 'error'

export const selectExperimentsErrorMessage = (state: RootState) => {
  const experiments = selectExperiments(state)
  if (experiments.status === 'error') {
    return experiments.message
  } else {
    throw new Error('experiments status is not error')
  }
}

export const selectExperimentList = (state: RootState) => {
  const experiments = selectExperiments(state)
  if (experiments.status === 'fulfilled') {
    return experiments.experimentList
  } else {
    throw new Error('experiments status is not fulfilled')
  }
}

export const selectExperimentUidList = (state: RootState) =>
  Object.keys(selectExperimentList(state))

export const selectExperiment = (uid: string) => (state: RootState) =>
  selectExperimentList(state)[uid]

export const selectExperimentTimeStamp = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).timestamp

export const selectExperimentName = (uid: string) => (state: RootState) =>
  selectExperiment(uid)(state).name

export const selectExperimentStatus =
  (uid: string) =>
  (state: RootState): EXPERIMENTS_STATUS => {
    const functions = selectExperimentList(state)[uid].functions
    const statusList = Object.values(functions).map((f) => f.status)
    if (statusList.findIndex((status) => status === 'error') >= 0) {
      return 'error'
    } else if (statusList.findIndex((status) => status === 'running') >= 0) {
      return 'running'
    } else {
      return 'success'
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

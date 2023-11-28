import {
  configureStore,
  ThunkAction,
  Action,
  combineReducers,
} from "@reduxjs/toolkit"

import {
  algorithmListReducer,
  algorithmNodeReducer,
  displayDataReducer,
  fileUploaderReducer,
  flowElementReducer,
  inputNodeReducer,
  nwbReducer,
  filesTreeReducer,
  handleTypeColorReducer,
  rightDrawerReducer,
  visualaizeItemReducer,
  snakemakeReducer,
  pipelineReducer,
  hdf5Reducer,
  experimentsReducer,
  workspaceReducer,
  userReducer,
} from "store/slice"

export const rootReducer = combineReducers({
  algorithmList: algorithmListReducer,
  algorithmNode: algorithmNodeReducer,
  displayData: displayDataReducer,
  fileUploader: fileUploaderReducer,
  flowElement: flowElementReducer,
  inputNode: inputNodeReducer,
  handleColor: handleTypeColorReducer,
  filesTree: filesTreeReducer,
  nwb: nwbReducer,
  rightDrawer: rightDrawerReducer,
  visualaizeItem: visualaizeItemReducer,
  snakemake: snakemakeReducer,
  pipeline: pipelineReducer,
  hdf5: hdf5Reducer,
  experiments: experimentsReducer,
  workspace: workspaceReducer,
  user: userReducer,
})

export const store = configureStore({
  reducer: rootReducer,
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware({
      serializableCheck: {
        ignoredActions: [
          "user/deleteUser/rejected",
          "user/updateUser/rejected",
          "workflow/fetchExperiment/rejected",
          "workflow/reproduceWorkflow/rejected",
          "workspace/getWorkspace/rejected",
          "workspace/getListUserShareWorkSpaces/rejected",
          "workspace/putWorkspaceList/rejected",
          "workspace/delWorkspaceList/rejected",
        ],
      },
    }),
})

export type AppDispatch = typeof store.dispatch
export type RootState = ReturnType<typeof store.getState>
export type AppThunk<ReturnType = void> = ThunkAction<
  ReturnType,
  RootState,
  unknown,
  Action<string>
>
export type ThunkApiConfig<T = unknown, PendingMeta = unknown> = {
  state: RootState
  dispatch: AppDispatch
  rejectValue: T
  penfingMeta: PendingMeta
}

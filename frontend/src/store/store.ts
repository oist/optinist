import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import { setupListeners } from '@reduxjs/toolkit/query'
import { webSocketApi } from '../api/Run/Run_old'
import {
  algorithmListReducer,
  algorithmNodeReducer,
  displayDataReducer,
  fileUploaderReducer,
  flowElementReducer,
  inputNodeReducer,
  runPipelineResultReducer,
  nwbReducer,
  filesTreeReducer,
  handleTypeColorReducer,
  rightDrawerReducer,
  visualaizeItemReducer,
  snakemakeReducer,
  pipelineReducer,
} from './slice'

export const store = configureStore({
  reducer: {
    algorithmList: algorithmListReducer,
    algorithmNode: algorithmNodeReducer,
    displayData: displayDataReducer,
    fileUploader: fileUploaderReducer,
    flowElement: flowElementReducer,
    inputNode: inputNodeReducer,
    handleColor: handleTypeColorReducer,
    filesTree: filesTreeReducer,
    nwb: nwbReducer,
    runPipelineResult: runPipelineResultReducer,
    rightDrawer: rightDrawerReducer,
    visualaizeItem: visualaizeItemReducer,
    snakemake: snakemakeReducer,
    pipeline: pipelineReducer,
    [webSocketApi.reducerPath]: webSocketApi.reducer,
  },
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware().concat(webSocketApi.middleware),
})

setupListeners(store.dispatch)

export type AppDispatch = typeof store.dispatch
export type RootState = ReturnType<typeof store.getState>
export type AppThunk<ReturnType = void> = ThunkAction<
  ReturnType,
  RootState,
  unknown,
  Action<string>
>
export type ThunkApiConfig<T = unknown> = {
  state: RootState
  dispatch: AppDispatch
  rejectValue: T
}

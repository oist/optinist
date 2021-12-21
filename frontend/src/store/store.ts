import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import { setupListeners } from '@reduxjs/toolkit/query'
import { webSocketApi } from '../api/Run/Run'
import {
  algorithmListReducer,
  algorithmNodeReducer,
  displayDataReducer,
  fileUploaderReducer,
  flowElementReducer,
  inputNodeReducer,
  layoutTabReducer,
  runPipelineResultReducer,
  nwbReducer,
  filesTreeReducer,
  handleTypeColorReducer,
} from './slice'

export const store = configureStore({
  reducer: {
    algorithmList: algorithmListReducer,
    algorithmNode: algorithmNodeReducer,
    displayData: displayDataReducer,
    fileUploader: fileUploaderReducer,
    flowElement: flowElementReducer,
    inputNode: inputNodeReducer,
    layoutTab: layoutTabReducer,
    handleColor: handleTypeColorReducer,
    filesTree: filesTreeReducer,
    nwb: nwbReducer,
    runPipelineResult: runPipelineResultReducer,
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

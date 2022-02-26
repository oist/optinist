import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import { setupListeners } from '@reduxjs/toolkit/query'
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
    rightDrawer: rightDrawerReducer,
    visualaizeItem: visualaizeItemReducer,
    snakemake: snakemakeReducer,
    pipeline: pipelineReducer,
    hdf5: hdf5Reducer,
    experiments: experimentsReducer,
  },
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
export type ThunkApiConfig<T = unknown, PendingMeta = unknown> = {
  state: RootState
  dispatch: AppDispatch
  rejectValue: T
  penfingMeta: PendingMeta
}

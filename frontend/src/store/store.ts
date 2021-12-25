import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import { setupListeners } from '@reduxjs/toolkit/query'
import { webSocketApi } from '../api/Run/Run'
import elementReducer from './slice/Element/Element'
import algorithmReducer from './slice/Algorithm/Algorithm'
import fileDataReducer from './slice/FileData/FileData'
import plotDataReducer from './slice/PlotData/PlotData'
import handleTypeColorReducer from './slice/HandleTypeColor/HandleTypeColor'
import filesTreeReducer from './slice/FilesTree/FilesTree'
import nwbReducer from './slice/NWB/NWB'

export const store = configureStore({
  reducer: {
    element: elementReducer,
    algorithm: algorithmReducer,
    fileData: fileDataReducer,
    plotData: plotDataReducer,
    handleColor: handleTypeColorReducer,
    filesTree: filesTreeReducer,
    nwb: nwbReducer,
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

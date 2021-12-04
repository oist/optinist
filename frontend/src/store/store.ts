import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import { setupListeners } from '@reduxjs/toolkit/query'
import elementReducer from './slice/Element/Element'
import algorithmReducer from './slice/Algorithm/Algorithm'
import uploadImageReducer from './slice/UploadImage/UploadImage'
import plotDataReducer from './slice/PlotData/PlotData'
import HandleTypeColorReducer from './slice/HandleTypeColor/HandleTypeColor'
import FilesReducer from './slice/Files/Files'
import { webSocketApi } from '../api/Run/Run'

export const store = configureStore({
  reducer: {
    element: elementReducer,
    algorithm: algorithmReducer,
    uploadImage: uploadImageReducer,
    plotData: plotDataReducer,
    handleColor: HandleTypeColorReducer,
    files: FilesReducer,
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

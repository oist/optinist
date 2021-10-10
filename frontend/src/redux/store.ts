import {
  configureStore,
  ThunkAction,
  Action,
  combineReducers,
} from '@reduxjs/toolkit'
import elementReducer from './slice/Element/Element'
import imageIndexReducer from './slice/ImageIndex/ImageIndex'
import outputReducer from './slice/Output/Output'

export const store = configureStore({
  reducer: combineReducers({
    element: elementReducer,
    imageIndex: imageIndexReducer,
    output: outputReducer,
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
export type ThunkApiConfig<T = unknown> = {
  state: RootState
  dispatch: AppDispatch
  rejectValue: T
}

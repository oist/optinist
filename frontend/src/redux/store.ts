import { configureStore, ThunkAction, Action } from '@reduxjs/toolkit'
import element from './slice/Element/Element'

export const store = configureStore({
  reducer: {
    element: element,
  },
})

export type AppDispatch = typeof store.dispatch
export type RootState = ReturnType<typeof store.getState>
export type AppThunk<ReturnType = void> = ThunkAction<
  ReturnType,
  RootState,
  unknown,
  Action<string>
>

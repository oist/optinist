import { createContext } from 'react'
import State from '../models/State'
import StateAction from '../models/StateAction'

const AppStateContext = createContext(
  {} as {
    state: State
    dispatch: React.Dispatch<StateAction>
  },
)

export default AppStateContext

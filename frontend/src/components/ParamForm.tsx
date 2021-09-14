import { useContext } from 'react'
import ParamItem from './ParamItem'

import AppStateContext from 'contexts/AppStateContext'

const ParamForm = () => {
  const { state } = useContext(AppStateContext)
  return (
    <div>
      <h2>{state.currentSelectedAlgo}</h2>
      {state.algorithms.map((algo) => {
        if (algo.name === state.currentSelectedAlgo) {
          return (
            <ul>
              {algo.parameters.map((param) => {
                return <ParamItem name={param.name} />
              })}
            </ul>
          )
        } else {
          return null
        }
      })}
      {/* <ParamItem name={'alpha'} />
      <ParamItem name={'beta'} />
      <ParamItem name={'gamma'} /> */}
    </div>
  )
}

export default ParamForm

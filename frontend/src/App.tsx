import { useReducer, useEffect } from 'react'
import './App.css'
import { Layout, Model, TabNode } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import { flexjson } from 'const/flexlayout'
import SideBar from 'components/SideBar'
import FlowChart from 'components/FlowChart'
import ParamForm from 'components/ParamForm'
import PlotOutput from 'components/PlotOutput'
import ImageViewer from 'components/ImageViewer'

import AppStateContext from 'contexts/AppStateContext'

import State from 'models/State'
import Algorithm from 'models/Algorithm'
import Parameter from 'models/Parameter'
import StateAction from 'models/StateAction'

const model = Model.fromJson(flexjson)

const reducer = (state: State, action: StateAction): State => {
  switch (action.type) {
    case 'AlgoSelect': {
      const updatedState: State = {
        currentSelectedAlgo: action.value,
        algorithms: state.algorithms,
      }
      return updatedState
    }
    case 'ParamUpdate': {
      const newAlgorithms: Algorithm[] = state.algorithms.map((algo) => {
        if (algo.name === state.currentSelectedAlgo) {
          const newParameters: Parameter[] = algo.parameters.map((param) => {
            if (param.name === action.param) {
              param.value = action.value
            }
            return param
          })
          algo.parameters = newParameters
        }
        return algo
      })

      const updatedState: State = {
        currentSelectedAlgo: state.currentSelectedAlgo,
        algorithms: newAlgorithms,
      }
      return updatedState
    }
    default:
      return state
  }
}

function App() {
  // アルゴリズム(CaImAn, suite2p...)ごとにパラメータの初期値を作っておく必要がある
  const initialParametersCaiman: Parameter[] = [
    { name: 'alpha_caiman', value: 30 },
    { name: 'beta_caiman', value: 30 },
    { name: 'gamma_caiman', value: 30 },
  ]
  const initialParametersS2P: Parameter[] = [
    { name: 'alpha_s2p', value: 30 },
    { name: 'beta_s2p', value: 30 },
    { name: 'gamma_s2p', value: 30 },
  ]
  // アルゴリズムと対応するパラメータの初期値を設定
  const initialAlgorithms: Algorithm[] = [
    { name: 'CaImAn', parameters: initialParametersCaiman },
    { name: 'Suite2P', parameters: initialParametersS2P },
    { name: 'algo3', parameters: initialParametersCaiman },
  ]

  // Appのstateの初期値
  const initialState: State = {
    currentSelectedAlgo: 'CaImAn',
    algorithms: initialAlgorithms,
  }

  const [state, dispatch] = useReducer(reducer, initialState)

  // 確認用
  useEffect(() => {
    console.log(state)
  }, [state])

  const factory = (node: TabNode) => {
    var component = node.getComponent()
    if (component === 'button') {
      return <button>{node.getName()}</button>
    } else if (component == 'flowchart') {
      return <FlowChart />
    } else if (component == 'sidebar') {
      return <SideBar />
    } else if (component == 'paramForm') {
      return <ParamForm />
    } else if (component == 'output') {
      return <PlotOutput />
    } else if (component == 'image') {
      return <ImageViewer />
    } else {
      return null
    }
  }

  return (
    <AppStateContext.Provider value={{ state, dispatch }}>
      <Layout model={model} factory={factory} />
    </AppStateContext.Provider>
  )
}

export default App

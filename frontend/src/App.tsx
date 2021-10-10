import './App.css'
import { Layout, Model, TabNode, Actions } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import { flexjson } from 'const/flexlayout'
import { SideBar } from 'components/TreeView'
import { FlowChart } from 'components/FlowChart'
import { ParamForm } from 'components/ParamForm/ParamForm'
import { PlotOutput } from 'components/PlotOutput'
import { ImageViewer } from 'components/ImageViewer'
import { ToolBar } from 'components/ToolBar'
import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import {
  runMassageSelector,
  runStatusSelector,
} from 'redux/slice/Element/ElementSelector'
import { RUN_STATUS } from 'redux/slice/Element/ElementType'

const model = Model.fromJson(flexjson)

function App() {
  const factory = (node: TabNode) => {
    var component = node.getComponent()
    if (component === 'button') {
      return <button>{node.getName()}</button>
    } else if (component === 'flowchart') {
      return <FlowChart />
    } else if (component === 'sidebar') {
      return <SideBar />
    } else if (component === 'paramForm') {
      return <ParamForm />
    } else if (component === 'output') {
      return <PlotOutput />
    } else if (component === 'image') {
      return <ImageViewer />
    } else {
      return null
    }
  }
  const runStatus = useSelector(runStatusSelector)
  React.useEffect(() => {
    if (runStatus === RUN_STATUS.SUCCESS) {
      model.doAction(Actions.selectTab('output0'))
    }
  }, [runStatus])
  return (
    <div id="container">
      <div className="app">
        <div className="toolbar">
          <ToolBar />
        </div>
        <div className="contents">
          <Layout model={model} factory={factory} />
        </div>
      </div>
    </div>
  )
}

export default App

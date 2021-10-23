import './App.css'
import { Layout, Model, TabNode } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import { flexjson } from 'const/flexlayout'
import { SideBar } from 'components/TreeView'
import { FlowChart } from 'components/FlowChart'
import { ParamForm } from 'components/ParamForm/ParamForm'
import { PlotOutput } from 'components/PlotOutput'
import { ImageViewer } from 'components/ImageViewer'
import { ToolBar } from 'components/ToolBar'
import React from 'react'

const model = Model.fromJson(flexjson)

export const FlexLayoutModelContext = React.createContext<Model>(model)
export const NodeIdContext = React.createContext<string>('')

function App() {
  const factory = (node: TabNode) => {
    var component = node.getComponent()
    const nodeId = node.getId().split('-')[0] // todo function化する
    switch (component) {
      case 'flowchart':
        return <FlowChart />
      case 'sidebar':
        return <SideBar />
      case 'paramForm':
        return (
          <NodeIdContext.Provider value={nodeId}>
            <ParamForm />
          </NodeIdContext.Provider>
        )
      case 'output':
        return (
          <NodeIdContext.Provider value={nodeId}>
            <PlotOutput />
          </NodeIdContext.Provider>
        )
      case 'image':
        return (
          <NodeIdContext.Provider value={nodeId}>
            <ImageViewer />
          </NodeIdContext.Provider>
        )
      default:
        return null
    }
  }

  return (
    <div id="container">
      <div className="app">
        <FlexLayoutModelContext.Provider value={model}>
          <div className="toolbar">
            <ToolBar />
          </div>
          <div className="contents">
            <Layout model={model} factory={factory} />
          </div>
        </FlexLayoutModelContext.Provider>
      </div>
    </div>
  )
}

export default App

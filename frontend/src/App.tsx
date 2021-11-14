import './App.css'
import { Layout, Model, TabNode } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import { flexjson } from 'const/flexlayout'
import { SideBar } from 'components/TreeView'
import { FlowChart } from 'components/FlowChart'
import { ParamForm } from 'components/ParamForm/ParamForm'
import { Plot } from 'components/Plot/Plot'
import { ToolBar } from 'components/ToolBar'
import React from 'react'
import { getNodeId, getSuffix } from 'utils/FlexLayoutUtils'
import { ImagePlot } from 'components/Plot/ImagePlot'

const model = Model.fromJson(flexjson)

export const FlexLayoutModelContext = React.createContext<Model>(model)
export const NodeIdContext = React.createContext<string>('')
export const OutputPlotContext = React.createContext<{
  nodeId: string
  outputKey: string
}>({ nodeId: '', outputKey: '' })

export const ImageDataContext = React.createContext<{
  nodeId: string
  outputKey: string | null // FileSelectorNodeでアップロードされた場合はoutputKeyがnull
}>({ nodeId: '', outputKey: null })

function App() {
  const factory = (node: TabNode) => {
    var component = node.getComponent()
    const layoutTabId = node.getId()
    const nodeId = getNodeId(layoutTabId)
    if (nodeId == null) {
      return null
    }
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
        const outputKey = getSuffix(layoutTabId)
        if (outputKey != null) {
          return (
            <OutputPlotContext.Provider value={{ nodeId, outputKey }}>
              <Plot />
            </OutputPlotContext.Provider>
          )
        } else {
          return null
        }
      case 'image':
        const key = getSuffix(layoutTabId)
        return (
          <ImageDataContext.Provider value={{ nodeId, outputKey: key }}>
            <ImagePlot />
          </ImageDataContext.Provider>
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

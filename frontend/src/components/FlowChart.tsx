import React, { useState, useContext, DragEvent } from 'react'
import 'style/flow.css'

import ReactFlow, {
  ReactFlowProvider,
  removeElements,
  addEdge,
  Controls,
  OnLoadParams,
  ElementId,
  Elements,
  Connection,
  Edge,
  Node,
} from 'react-flow-renderer'
import { initialElements } from 'const/flowchart'
import AppStateContext from 'contexts/AppStateContext'

const onDragOver = (event: DragEvent) => {
  event.preventDefault()
  event.dataTransfer.dropEffect = 'move'
}

let id = 0
const getId = (): ElementId => `dndnode_${id++}`

const FlowChart = () => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const [elements, setElements] = useState<Elements>(initialElements)
  const { dispatch } = useContext(AppStateContext)

  const onConnect = (params: Connection | Edge) =>
    setElements((els) =>
      addEdge({ ...params, type: 'smoothstep', animated: false }, els),
    )

  const onElementClick = (
    event: React.MouseEvent<Element, MouseEvent>,
    element: any,
  ) => {
    if (event.isTrusted) {
      dispatch({ type: 'AlgoSelect', value: element.data.label })
    }
  }

  const onElementsRemove = (elementsToRemove: Elements) =>
    setElements((els) => removeElements(elementsToRemove, els))

  const onLoad = (_reactFlowInstance: OnLoadParams) =>
    setReactFlowInstance(_reactFlowInstance)

  const onDrop = (event: DragEvent) => {
    event.preventDefault()

    if (reactFlowInstance) {
      const name = event.dataTransfer.getData('application/reactflow')
      const position = reactFlowInstance.project({
        x: event.clientX - 50 - 250,
        y: event.clientY - 50,
      })

      var type = 'default'
      if (name.includes('data')) {
        type = 'input'
      } else if (name.includes('output')) {
        type = 'output'
      }

      const newNode: Node = {
        id: getId(),
        type: type,
        position,
        data: { label: `${name}` },
      }

      setElements((es) => es.concat(newNode))
    }
  }

  return (
    <div className="flow">
      <ReactFlowProvider>
        <div className="reactflow-wrapper">
          <ReactFlow
            elements={elements}
            onElementClick={onElementClick}
            onElementsRemove={onElementsRemove}
            onConnect={onConnect}
            onLoad={onLoad}
            onDrop={onDrop}
            onDragOver={onDragOver}
          >
            <Controls />
          </ReactFlow>
        </div>
      </ReactFlowProvider>
    </div>
  )
}

export default FlowChart

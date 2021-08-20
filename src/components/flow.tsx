import { useState, DragEvent } from 'react'
import Siderbar from './Siderbar'
import './dnd.css'

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

const initialElements = [
  {
    id: '1',
    type: 'input',
    data: { label: 'input node' },
    position: { x: 250, y: 5 },
  },
]

const onDragOver = (event: DragEvent) => {
  event.preventDefault()
  event.dataTransfer.dropEffect = 'move'
}

let id = 0
const getId = (): ElementId => `dndnode_${id++}`

const BasicFlow = () => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const [elements, setElements] = useState<Elements>(initialElements)

  const onConnect = (params: Connection | Edge) =>
    setElements((els) => addEdge(params, els))
  const onElementsRemove = (elementsToRemove: Elements) =>
    setElements((els) => removeElements(elementsToRemove, els))
  const onLoad = (_reactFlowInstance: OnLoadParams) =>
    setReactFlowInstance(_reactFlowInstance)

  const onDrop = (event: DragEvent) => {
    event.preventDefault()

    if (reactFlowInstance) {
      const type = event.dataTransfer.getData('application/reactflow')
      const position = reactFlowInstance.project({
        x: event.clientX,
        y: event.clientY - 40,
      })
      const newNode: Node = {
        id: getId(),
        type,
        position,
        data: { label: `${type} node` },
      }

      setElements((es) => es.concat(newNode))
    }
  }

  return (
    <div className="dndflow">
      <Siderbar />
      <ReactFlowProvider>
        <div className="reactflow-wrapper">
          <ReactFlow
            elements={elements}
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

export default BasicFlow

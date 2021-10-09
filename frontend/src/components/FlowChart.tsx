import React, { useState, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
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
import { setFlowElements, addFlowElement } from 'redux/slice/Element/Element'
import { flowElementsSelector } from 'redux/slice/Element/ElementSelector'
import { NodeDataType, NodeType } from 'redux/slice/Element/ElementType'
import 'style/flow.css'
import { FileSelectorNode } from './FileSelectorNode'
import { clickNode } from 'redux/slice/Element/ElementAction'
import { isNodeData } from 'redux/slice/Element/ElementUtils'

let id = 0
const getId = (): ElementId => `dndnode_${id++}`

export const FlowChart = React.memo(() => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const flowElements = useSelector(flowElementsSelector)
  const dispatch = useDispatch()

  const nodeTypes = {
    selectorNode: FileSelectorNode,
  }

  const onConnect = (params: Connection | Edge) => {
    dispatch(
      setFlowElements(
        addEdge(
          { ...params, type: 'smoothstep', animated: false },
          flowElements,
        ),
      ),
    )
  }

  const onElementClick = (
    event: React.MouseEvent<Element, MouseEvent>,
    element: Node<NodeDataType> | Edge<any>,
  ) => {
    if (event.isTrusted && isNodeData(element) && element.data) {
      dispatch(clickNode({ id: element.id, type: element.data.type }))
    }
  }

  const onElementsRemove = (elementsToRemove: Elements) => {
    dispatch(setFlowElements(removeElements(elementsToRemove, flowElements)))
  }

  const onLoad = (_reactFlowInstance: OnLoadParams) =>
    setReactFlowInstance(_reactFlowInstance)

  const onDragOver = (event: DragEvent) => {
    event.preventDefault()
    event.dataTransfer.dropEffect = 'move'
  }

  const onDrop = (event: DragEvent) => {
    event.preventDefault()

    if (reactFlowInstance) {
      const name = event.dataTransfer.getData('application/reactflow')
      const position = reactFlowInstance.project({
        x: event.clientX - 50 - 250,
        y: event.clientY - 50,
      })

      let type: NodeType = 'algo'
      if (name.includes('data')) {
        type = 'input'
      } else if (name.includes('output')) {
        type = 'output'
      }

      const newNode: Node<NodeDataType> = {
        id: getId(),
        type: type,
        position,
        data: { label: name, type },
      }

      dispatch(addFlowElement(newNode))
    }
  }

  return (
    <div className="flow">
      <ReactFlowProvider>
        <div className="reactflow-wrapper">
          <ReactFlow
            elements={flowElements}
            onElementClick={onElementClick}
            onElementsRemove={onElementsRemove}
            onConnect={onConnect}
            onLoad={onLoad}
            onDrop={onDrop}
            onDragOver={onDragOver}
            nodeTypes={nodeTypes}
          >
            <Controls />
          </ReactFlow>
        </div>
      </ReactFlowProvider>
    </div>
  )
})

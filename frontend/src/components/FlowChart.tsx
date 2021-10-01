import React, { useState, DragEvent } from 'react'
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

import ColorSelectorNode from './FileSelectorNode'
import {
  flowElementsSelector,
  algoParamsSelector,
} from 'redux/slice/Element/ElementSelector'
import { useSelector, useDispatch } from 'react-redux'
import { setFlowElements, setCurrentElement } from 'redux/slice/Element/Element'

let id = 0
const getId = (): ElementId => `dndnode_${id++}`

export const FlowChart = React.memo(() => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const flowElements = useSelector(flowElementsSelector)
  const algoParams = useSelector(algoParamsSelector)
  const dispatch = useDispatch()

  const nodeTypes = {
    selectorNode: ColorSelectorNode,
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
    element: any,
  ) => {
    if (event.isTrusted) {
      dispatch(setCurrentElement(element.data.label))
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

      dispatch(setFlowElements(flowElements.concat(newNode)))
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

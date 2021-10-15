import React, { useState, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import ReactFlow, {
  ReactFlowProvider,
  removeElements,
  addEdge,
  Controls,
  OnLoadParams,
  Elements,
  Connection,
  Edge,
  Node,
} from 'react-flow-renderer'
import { setFlowElements, addFlowElement } from 'redux/slice/Element/Element'
import {
  flowElementsSelector,
  maxElementIdSelector,
} from 'redux/slice/Element/ElementSelector'
import { NodeData, NODE_DATA_TYPE, NODE_DATA_TYPE_SET } from 'const/NodeData'
import 'style/flow.css'
import { FileSelectorNode } from './FileSelectorNode'
import { clickNode } from 'redux/slice/Element/ElementAction'
import { isNodeData } from 'utils/ElementUtils'

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
    element: Node<NodeData> | Edge<any>,
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

  const maxElementId = useSelector(maxElementIdSelector)
  const onDrop = (event: DragEvent) => {
    event.preventDefault()

    if (reactFlowInstance) {
      const name = event.dataTransfer.getData('application/reactflow')
      const position = reactFlowInstance.project({
        x: event.clientX - 50 - 250,
        y: event.clientY - 50,
      })

      let nodeType = 'default'
      let dataType: NODE_DATA_TYPE = 'algo'
      if (name.includes('data')) {
        dataType = NODE_DATA_TYPE_SET.DATA
        nodeType = 'input'
      } else if (name.includes('output')) {
        dataType = NODE_DATA_TYPE_SET.OUTPUT
        nodeType = 'output'
      }

      const newNode: Node<NodeData> = {
        id: String(maxElementId + 1),
        type: nodeType,
        position,
        data: { label: name, type: dataType },
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

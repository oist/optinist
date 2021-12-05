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
import { setFlowElements, addFlowElement } from 'store/slice/Element/Element'
import {
  flowElementsSelector,
  maxElementIdSelector,
} from 'store/slice/Element/ElementSelector'
import { NodeData, NODE_DATA_TYPE, NODE_DATA_TYPE_SET } from 'const/NodeData'
import 'style/flow.css'
import { ImageFileNode } from './FlowChartNode/ImageFileNode'
import { AlgorithmNode } from './FlowChartNode/AlgorithmNode'
import { CsvFileNode } from './FlowChartNode/CsvFileNode'
import { FlexLayoutModelContext } from 'App'
import { useGetDeleteTabActions } from 'FlexLayoutHook'

export const FlowChart = React.memo(() => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const flowElements = useSelector(flowElementsSelector)
  const dispatch = useDispatch()

  const nodeTypes = {
    ImageFileNode,
    CsvFileNode,
    AlgorithmNode,
  }

  const onConnect = (params: Connection | Edge) => {
    dispatch(
      setFlowElements(
        addEdge(
          {
            ...params,
            type: 'smoothstep',
            animated: false,
            style: { width: 5 },
          },
          flowElements,
        ),
      ),
    )
  }

  const model = React.useContext(FlexLayoutModelContext)
  const getDeleteTabActions = useGetDeleteTabActions()
  const onElementsRemove = (elementsToRemove: Elements) => {
    dispatch(setFlowElements(removeElements(elementsToRemove, flowElements)))
    const ids = elementsToRemove.map((elm) => elm.id)
    const actions = getDeleteTabActions(...ids)
    actions.forEach((action) => {
      model.doAction(action)
    })
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
      const type = event.dataTransfer.getData('type')
      const name = event.dataTransfer.getData('nodeName')
      const path = event.dataTransfer.getData('path')
      const position = reactFlowInstance.project({
        x: event.clientX - 50 - 250,
        y: event.clientY - 100,
      })
      let nodeType = ''
      let dataType: NODE_DATA_TYPE = 'algo'
      switch (type) {
        case NODE_DATA_TYPE_SET.ALGO:
          dataType = NODE_DATA_TYPE_SET.ALGO
          nodeType = 'AlgorithmNode'
          break
        case NODE_DATA_TYPE_SET.IMAGE:
          dataType = NODE_DATA_TYPE_SET.IMAGE
          nodeType = 'ImageFileNode'
          break
        case NODE_DATA_TYPE_SET.CSV:
          dataType = NODE_DATA_TYPE_SET.CSV
          nodeType = 'CsvFileNode'
          break
        default:
          break
      }

      const newNode: Node<NodeData> = {
        id: String(maxElementId + 1),
        type: nodeType,
        position,
        data: { label: name, type: dataType, path },
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
            // onElementClick={onElementClick}
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

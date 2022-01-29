import React, { useState, DragEvent, MouseEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import ReactFlow, {
  ReactFlowProvider,
  addEdge,
  Controls,
  OnLoadParams,
  Elements,
  Connection,
  Edge,
  Node,
} from 'react-flow-renderer'

import 'style/flow.css'
import {
  NodeData,
  NODE_TYPE_SET,
} from 'store/slice/FlowElement/FlowElementType'
import { FILE_TYPE, FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  addFlowElementNode,
  deleteFlowElements,
  editFlowElementPositionById,
  setFlowElements,
} from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectFlowElements,
  selectMaxElementId,
} from 'store/slice/FlowElement/FlowElementSelectors'
import { ImageFileNode } from './FlowChartNode/ImageFileNode'
import { AlgorithmNode } from './FlowChartNode/AlgorithmNode'
import { CsvFileNode } from './FlowChartNode/CsvFileNode'
import { NWBFileNode } from './FlowChartNode/NWBFileNode'
import { CustomEdge } from './CustomEdge'

const componentTypes = {
  ImageFileNode,
  CsvFileNode,
  NWBFileNode,
  AlgorithmNode,
} as const

const edgeTypes = {
  buttonedge: CustomEdge,
} as const

export const ReactFlowComponent = React.memo(() => {
  const [reactFlowInstance, setReactFlowInstance] = useState<OnLoadParams>()
  const flowElements = useSelector(selectFlowElements)
  const dispatch = useDispatch()

  const onConnect = (params: Connection | Edge) => {
    dispatch(
      setFlowElements(
        addEdge(
          {
            ...params,
            animated: false,
            style: { width: 5 },
            type: 'buttonedge',
          },
          flowElements,
        ),
      ),
    )
  }

  const onElementsRemove = (elementsToRemove: Elements) => {
    dispatch(deleteFlowElements(elementsToRemove))
  }

  const onLoad = (_reactFlowInstance: OnLoadParams) =>
    setReactFlowInstance(_reactFlowInstance)

  const onDragOver = (event: DragEvent) => {
    event.preventDefault()
    event.dataTransfer.dropEffect = 'move'
  }

  const maxElementId = useSelector(selectMaxElementId)
  const onDrop = (event: DragEvent) => {
    event.preventDefault()

    if (reactFlowInstance) {
      const position = reactFlowInstance.project({
        x: event.clientX - 50 - 250,
        y: event.clientY - 100,
      })
      const name = event.dataTransfer.getData('nodeName')
      const nodeType = event.dataTransfer.getData('nodeType')
      if (nodeType === NODE_TYPE_SET.INPUT) {
        let fileType: FILE_TYPE = FILE_TYPE_SET.CSV
        let componentType = ''
        switch (event.dataTransfer.getData('fileType')) {
          case FILE_TYPE_SET.CSV:
            componentType = 'CsvFileNode'
            break
          case FILE_TYPE_SET.IMAGE:
            componentType = 'ImageFileNode'
            fileType = FILE_TYPE_SET.IMAGE
            break
          case FILE_TYPE_SET.NWB:
            componentType = 'NWBFileNode'
            fileType = FILE_TYPE_SET.NWB
            break
        }

        const newNode: Node<NodeData> = {
          id: String(maxElementId + 1),
          type: componentType,
          position,
          data: { label: name, type: nodeType },
        }
        dispatch(
          addFlowElementNode({ node: newNode, inputNodeInfo: { fileType } }),
        )
      } else if (nodeType === NODE_TYPE_SET.ALGORITHM) {
        const functionPath = event.dataTransfer.getData('functionPath')
        const newNode: Node<NodeData> = {
          id: String(maxElementId + 1),
          type: 'AlgorithmNode',
          position,
          data: { label: name, type: nodeType },
        }
        dispatch(
          addFlowElementNode({
            node: newNode,
            algoNodeInfo: { name, functionPath },
          }),
        )
      }
    }
  }

  const onMoveEnd = (event: MouseEvent, node: Node) => {
    dispatch(
      editFlowElementPositionById({
        nodeId: node.id,
        coord: { x: node.position.x, y: node.position.y },
      }),
    )
  }

  return (
    <div className="flow">
      <ReactFlowProvider>
        <div className="reactflow-wrapper">
          <ReactFlow
            elements={flowElements}
            onElementsRemove={onElementsRemove}
            onConnect={onConnect}
            onLoad={onLoad}
            onDrop={onDrop}
            onDragOver={onDragOver}
            onNodeDragStop={onMoveEnd}
            nodeTypes={componentTypes}
            edgeTypes={edgeTypes}
          >
            <Controls />
          </ReactFlow>
        </div>
      </ReactFlowProvider>
    </div>
  )
})

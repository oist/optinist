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
  FlowTransform,
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
  setFlowPosition,
} from 'store/slice/FlowElement/FlowElementSlice'
import {
  selectFlowElements,
  selectFlowPosition,
} from 'store/slice/FlowElement/FlowElementSelectors'
import { ImageFileNode } from './FlowChartNode/ImageFileNode'
import { CsvFileNode } from './FlowChartNode/CsvFileNode'
import { HDF5FileNode } from './FlowChartNode/HDF5FileNode'
import { AlgorithmNode } from './FlowChartNode/AlgorithmNode'
import { CustomEdge } from './CustomEdge'
import { nanoid } from '@reduxjs/toolkit'

const componentTypes = {
  ImageFileNode,
  CsvFileNode,
  HDF5FileNode,
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

  const onDrop = (event: DragEvent) => {
    event.preventDefault()
    const id = nanoid()
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
            fileType = FILE_TYPE_SET.CSV
            break
          case FILE_TYPE_SET.IMAGE:
            componentType = 'ImageFileNode'
            fileType = FILE_TYPE_SET.IMAGE
            break
          case FILE_TYPE_SET.HDF5:
            componentType = 'HDF5FileNode'
            fileType = FILE_TYPE_SET.HDF5
            break
        }

        const newNode: Node<NodeData> = {
          id,
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
          id,
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

  const onNodeDragStop = (event: MouseEvent, node: Node) => {
    dispatch(
      editFlowElementPositionById({
        nodeId: node.id,
        coord: { x: node.position.x, y: node.position.y },
      }),
    )
  }

  const flowPosition = useSelector(selectFlowPosition)

  const onMoveEnd = (event: FlowTransform | undefined) => {
    if (event !== undefined) {
      dispatch(setFlowPosition(event))
    }
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
            onNodeDragStop={onNodeDragStop}
            nodeTypes={componentTypes}
            edgeTypes={edgeTypes}
            defaultPosition={[flowPosition.x, flowPosition.y]}
            defaultZoom={flowPosition.zoom}
            onMoveEnd={onMoveEnd}
          >
            <Controls />
          </ReactFlow>
        </div>
      </ReactFlowProvider>
    </div>
  )
})

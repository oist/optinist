import React, { useEffect, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { nanoid } from '@reduxjs/toolkit'
import { TreeView, TreeItem } from '@material-ui/lab'
import { makeStyles, Typography } from '@material-ui/core'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import IconButton from '@material-ui/core/IconButton'
import AddIcon from '@material-ui/icons/Add'

import {
  selectAlgorithmListIsLated,
  selectAlgorithmListTree,
} from 'store/slice/AlgorithmList/AlgorithmListSelectors'
import { AlgorithmNodeType } from 'store/slice/AlgorithmList/AlgorithmListType'
import { isAlgoChild } from 'store/slice/AlgorithmList/AlgorithmListUtils'
import { getAlgoList } from 'store/slice/AlgorithmList/AlgorithmListActions'
import { FILE_TYPE, FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  NodeData,
  NODE_TYPE,
  NODE_TYPE_SET,
} from 'store/slice/FlowElement/FlowElementType'
import {
  addFlowElementNode,
  setElementCoord,
} from 'store/slice/FlowElement/FlowElementSlice'
import { Node } from 'react-flow-renderer'
import { selectElementCoord } from 'store/slice/FlowElement/FlowElementSelectors'

const useStyles = makeStyles({
  root: {
    flexGrow: 1,
    height: '100%',
  },
})

export const AlgorithmTreeView = React.memo(() => {
  const dispatch = useDispatch()
  const classes = useStyles()
  const algoList = useSelector(selectAlgorithmListTree)
  const isLatest = useSelector(selectAlgorithmListIsLated)

  useEffect(() => {
    if (!isLatest) {
      dispatch(getAlgoList())
    }
  }, [dispatch, isLatest])

  return (
    <TreeView
      className={classes.root}
      defaultCollapseIcon={<ExpandMoreIcon />}
      defaultExpandIcon={<ChevronRightIcon />}
    >
      <TreeItem nodeId="Data" label="Data">
        <InputNodeComponent
          fileName={'image'}
          nodeName={'imageData'}
          fileType={FILE_TYPE_SET.IMAGE}
        />
        <InputNodeComponent
          fileName={'csv'}
          nodeName={'csvData'}
          fileType={FILE_TYPE_SET.CSV}
        />
        <InputNodeComponent
          fileName={'nwb'}
          nodeName={'nwbData'}
          fileType={FILE_TYPE_SET.NWB}
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponent name={name} node={node} key={i.toFixed()} />
        ))}
      </TreeItem>
    </TreeView>
  )
})

const InputNodeComponent = React.memo<{
  fileName: string
  nodeName: string
  fileType: FILE_TYPE
}>(({ fileName, nodeName, fileType }) => {
  const dispatch = useDispatch()
  const elementCoord = useSelector(selectElementCoord)

  const onDataNodeClick = (
    nodeType: NODE_TYPE,
    nodeName: string,
    fileType: FILE_TYPE,
  ) => {
    const position = {
      x: elementCoord.x,
      y: elementCoord.y,
    }
    dispatch(setElementCoord())

    let componentType = ''
    switch (fileType) {
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
      id: nanoid(),
      type: componentType,
      position,
      data: { label: nodeName, type: nodeType },
    }
    dispatch(addFlowElementNode({ node: newNode, inputNodeInfo: { fileType } }))
  }

  const onDataNodeDragStart = (
    event: DragEvent,
    nodeName: string,
    fileType: FILE_TYPE,
  ) => {
    if (event.dataTransfer != null) {
      event.dataTransfer.setData('nodeName', nodeName)
      event.dataTransfer.setData('nodeType', NODE_TYPE_SET.INPUT)
      event.dataTransfer.setData('fileType', fileType)
      event.dataTransfer.effectAllowed = 'move'
    }
  }

  return (
    <TreeItem
      nodeId={fileName}
      label={
        <AddButton
          name={fileName}
          onClick={() =>
            onDataNodeClick(NODE_TYPE_SET.INPUT, nodeName, fileType)
          }
        />
      }
      onDragStart={(event: DragEvent) =>
        onDataNodeDragStart(event, nodeName, fileType)
      }
      draggable
    ></TreeItem>
  )
})

const AlgoNodeComponent = React.memo<{
  name: string
  node: AlgorithmNodeType
}>(({ name, node }) => {
  const dispatch = useDispatch()
  const elementCoord = useSelector(selectElementCoord)

  const onAlgoNodeClick = (nodeName: string, functionPath: string) => {
    const name = nodeName
    const position = {
      x: elementCoord.x,
      y: elementCoord.y,
    }
    dispatch(setElementCoord())

    const newNode: Node<NodeData> = {
      id: nanoid(),
      type: 'AlgorithmNode',
      position,
      data: { label: name, type: NODE_TYPE_SET.ALGORITHM },
    }
    dispatch(
      addFlowElementNode({
        node: newNode,
        algoNodeInfo: { name, functionPath },
      }),
    )
  }

  const onAlgoNodeDragStart = (
    event: DragEvent,
    nodeName: string,
    functionPath: string,
  ) => {
    if (event.dataTransfer != null) {
      event.dataTransfer.setData('nodeName', nodeName)
      event.dataTransfer.setData('nodeType', NODE_TYPE_SET.ALGORITHM)
      event.dataTransfer.setData('functionPath', functionPath)
      event.dataTransfer.effectAllowed = 'move'
    }
  }

  if (isAlgoChild(node)) {
    return (
      <TreeItem
        nodeId={name}
        label={
          <AddButton
            name={name}
            onClick={() => onAlgoNodeClick(name, node.functionPath)}
          />
        }
        onDragStart={(event: DragEvent) =>
          onAlgoNodeDragStart(event, name, node.functionPath)
        }
        draggable
      />
    )
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        {Object.entries(node.children).map(([name, node], i) => (
          <AlgoNodeComponent name={name} node={node} key={i.toFixed()} />
        ))}
      </TreeItem>
    )
  }
})

const AddButton = React.memo<{
  name: string
  onClick: () => void
}>(({ name, onClick }) => {
  return (
    <>
      <IconButton aria-label="add" style={{ padding: 2 }}>
        <AddIcon onClick={() => onClick()} />
      </IconButton>
      <Typography
        variant="inherit"
        style={{
          textOverflow: 'ellipsis',
          overflow: 'visible',
          width: '8rem',
          display: 'inline-block',
        }}
      >
        {name}
      </Typography>
    </>
  )
})

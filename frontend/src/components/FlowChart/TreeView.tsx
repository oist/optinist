import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { nanoid } from '@reduxjs/toolkit'
import TreeView from '@mui/lab/TreeView'
import TreeItem, { treeItemClasses } from '@mui/lab/TreeItem'
import { styled, Typography } from '@mui/material'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'
import ChevronRightIcon from '@mui/icons-material/ChevronRight'
import IconButton from '@mui/material/IconButton'
import AddIcon from '@mui/icons-material/Add'

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

export const AlgorithmTreeView = React.memo(() => {
  const dispatch = useDispatch()
  const algoList = useSelector(selectAlgorithmListTree)
  const isLatest = useSelector(selectAlgorithmListIsLated)

  useEffect(() => {
    if (!isLatest) {
      dispatch(getAlgoList())
    }
  }, [dispatch, isLatest])

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

  return (
    <TreeView
      sx={{
        flexGrow: 1,
        height: '100%',
      }}
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
          fileName={'hdf5'}
          nodeName={'hdf5Data'}
          fileType={FILE_TYPE_SET.HDF5}
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponent
            name={name}
            node={node}
            onAlgoNodeClick={(name, functionPath) =>
              onAlgoNodeClick(name, functionPath)
            }
            key={i.toFixed()}
          />
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
      case FILE_TYPE_SET.HDF5:
        componentType = 'HDF5FileNode'
        fileType = FILE_TYPE_SET.HDF5
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

  return (
    <LeafItem
      nodeId={fileName}
      label={
        <AddButton
          name={fileName}
          onClick={() =>
            onDataNodeClick(NODE_TYPE_SET.INPUT, nodeName, fileType)
          }
        />
      }
    />
  )
})

const AlgoNodeComponent = React.memo<{
  name: string
  node: AlgorithmNodeType
  onAlgoNodeClick: (nodeName: string, functionPath: string) => void
}>(({ name, node, onAlgoNodeClick }) => {
  if (isAlgoChild(node)) {
    return (
      <LeafItem
        nodeId={name}
        label={
          <AddButton
            name={name}
            onClick={() => onAlgoNodeClick(name, node.functionPath)}
          />
        }
      />
    )
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        {Object.entries(node.children).map(([name, node], i) => (
          <AlgoNodeComponent
            name={name}
            node={node}
            key={i.toFixed()}
            onAlgoNodeClick={(name, functionPath) =>
              onAlgoNodeClick(name, functionPath)
            }
          />
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
      <IconButton aria-label="add" style={{ padding: 2 }} size="large">
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

// 未使用icon分の幅を消す
const LeafItem = styled(TreeItem)({
  [`& .${treeItemClasses.iconContainer}`]: {
    margin: 0,
    width: 0,
  },
})

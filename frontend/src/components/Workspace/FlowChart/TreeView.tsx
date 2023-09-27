import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import TreeView from '@mui/lab/TreeView'
import TreeItem, { treeItemClasses } from '@mui/lab/TreeItem'
import { styled, Typography } from '@mui/material'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'
import ChevronRightIcon from '@mui/icons-material/ChevronRight'
import IconButton from '@mui/material/IconButton'
import AddIcon from '@mui/icons-material/Add'
import { useDrag } from 'react-dnd'

import {
  selectAlgorithmListIsLated,
  selectAlgorithmListTree,
} from 'store/slice/AlgorithmList/AlgorithmListSelectors'
import {
  AlgorithmChild,
  AlgorithmNodeType,
} from 'store/slice/AlgorithmList/AlgorithmListType'
import { isAlgoChild } from 'store/slice/AlgorithmList/AlgorithmListUtils'
import { getAlgoList } from 'store/slice/AlgorithmList/AlgorithmListActions'
import { FILE_TYPE, FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'
import {
  NODE_TYPE,
  NODE_TYPE_SET,
} from 'store/slice/FlowElement/FlowElementType'
import {
  addAlgorithmNode,
  addInputNode,
} from 'store/slice/FlowElement/FlowElementActions'
import { getNanoId } from 'utils/nanoid/NanoIdUtils'
import { REACT_FLOW_NODE_TYPE, REACT_FLOW_NODE_TYPE_KEY } from 'const/flowchart'
import {
  DND_ITEM_TYPE_SET,
  TreeItemCollectedProps,
  TreeItemDragObject,
  TreeItemDropResult,
} from './DnDItemType'

export const AlgorithmTreeView = React.memo(() => {
  const dispatch = useDispatch()
  const algoList = useSelector(selectAlgorithmListTree)
  const isLatest = useSelector(selectAlgorithmListIsLated)

  useEffect(() => {
    if (!isLatest) {
      dispatch(getAlgoList())
    }
  }, [dispatch, isLatest])

  const onAddAlgoNode = React.useCallback(
    (
      nodeName: string,
      functionPath: string,
      position?: { x: number; y: number },
    ) => {
      const name = nodeName
      const newNode = {
        id: `${name}_${getNanoId()}`,
        type: REACT_FLOW_NODE_TYPE_KEY.AlgorithmNode,
        data: { label: name, type: NODE_TYPE_SET.ALGORITHM },
        position,
      }
      dispatch(
        addAlgorithmNode({
          node: newNode,
          name,
          functionPath,
        }),
      )
    },
    [dispatch],
  )

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
        <InputNodeComponent
          fileName={'fluo'}
          nodeName={'fluoData'}
          fileType={FILE_TYPE_SET.FLUO}
        />
        <InputNodeComponent
          fileName={'behavior'}
          nodeName={'behaviorData'}
          fileType={FILE_TYPE_SET.BEHAVIOR}
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponentRecursive
            name={name}
            node={node}
            onAddAlgoNode={onAddAlgoNode}
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

  const onAddDataNode = React.useCallback(
    (
      nodeType: NODE_TYPE,
      nodeName: string,
      fileType: FILE_TYPE,
      position?: { x: number; y: number },
    ) => {
      let reactFlowNodeType: REACT_FLOW_NODE_TYPE | '' = ''
      switch (fileType) {
        case FILE_TYPE_SET.CSV:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.CsvFileNode
          break
        case FILE_TYPE_SET.IMAGE:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.ImageFileNode
          fileType = FILE_TYPE_SET.IMAGE
          break
        case FILE_TYPE_SET.HDF5:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.HDF5FileNode
          fileType = FILE_TYPE_SET.HDF5
          break
        case FILE_TYPE_SET.FLUO:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.FluoFileNode
          fileType = FILE_TYPE_SET.FLUO
          break
        case FILE_TYPE_SET.BEHAVIOR:
          reactFlowNodeType = REACT_FLOW_NODE_TYPE_KEY.BehaviorFileNode
          fileType = FILE_TYPE_SET.BEHAVIOR
          break
      }
      const newNode = {
        id: `input_${getNanoId()}`,
        type: reactFlowNodeType,
        data: { label: nodeName, type: nodeType },
        position,
      }
      dispatch(addInputNode({ node: newNode, fileType }))
    },
    [dispatch],
  )

  const { isDragging, dragRef } = useLeafItemDrag(
    React.useCallback(
      (position) => {
        onAddDataNode(NODE_TYPE_SET.INPUT, nodeName, fileType, position)
      },
      [onAddDataNode, nodeName, fileType],
    ),
  )

  return (
    <LeafItem
      ref={dragRef}
      style={{
        opacity: isDragging ? 0.6 : 1,
      }}
      onFocusCapture={(e) => e.stopPropagation()}
      nodeId={fileName}
      label={
        <AddButton
          name={fileName}
          onClick={() => onAddDataNode(NODE_TYPE_SET.INPUT, nodeName, fileType)}
        />
      }
    />
  )
})

const AlgoNodeComponentRecursive = React.memo<{
  name: string
  node: AlgorithmNodeType
  onAddAlgoNode: (
    nodeName: string,
    functionPath: string,
    position?: { x: number; y: number },
  ) => void
}>(({ name, node, onAddAlgoNode }) => {
  if (isAlgoChild(node)) {
    return (
      <AlgoNodeComponent
        name={name}
        node={node}
        onAddAlgoNode={onAddAlgoNode}
      />
    )
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        {Object.entries(node.children).map(([name, node], i) => (
          <AlgoNodeComponentRecursive
            name={name}
            node={node}
            key={i.toFixed()}
            onAddAlgoNode={onAddAlgoNode}
          />
        ))}
      </TreeItem>
    )
  }
})

const AlgoNodeComponent = React.memo<{
  name: string
  node: AlgorithmChild
  onAddAlgoNode: (
    nodeName: string,
    functionPath: string,
    position?: { x: number; y: number },
  ) => void
}>(({ name, node, onAddAlgoNode }) => {
  const { isDragging, dragRef } = useLeafItemDrag(
    React.useCallback(
      (position) => {
        onAddAlgoNode(name, node.functionPath, position)
      },
      [onAddAlgoNode, name, node],
    ),
  )
  return (
    <LeafItem
      ref={dragRef}
      style={{
        opacity: isDragging ? 0.6 : 1,
      }}
      onFocusCapture={(e) => e.stopPropagation()}
      nodeId={name}
      label={
        <AddButton
          name={name}
          onClick={() => onAddAlgoNode(name, node.functionPath)}
        />
      }
    />
  )
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
  // background: 'red',
  [`& .${treeItemClasses.iconContainer}`]: {
    margin: 0,
    width: 0,
  },
})

function useLeafItemDrag(
  onDragEnd: (positon: { x: number; y: number }) => void,
) {
  const [{ isDragging }, dragRef] = useDrag<
    TreeItemDragObject,
    TreeItemDropResult,
    TreeItemCollectedProps
  >(
    () => ({
      type: DND_ITEM_TYPE_SET.TREE_ITEM,
      end: (_, monitor) => {
        const position = monitor.getDropResult()?.position
        if (monitor.didDrop() && position != null) {
          onDragEnd(position)
        }
      },
      collect: (monitor) => ({
        isDragging: monitor.isDragging(),
      }),
    }),
    [onDragEnd],
  )
  return { isDragging, dragRef }
}

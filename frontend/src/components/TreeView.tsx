import React, { useEffect, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

import { algoListSelector } from 'store/slice/Algorithm/AlgorithmSelector'
import { AlgoListType, AlgoNodeType } from 'store/slice/Algorithm/AlgorithmType'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { getAlgoList } from 'store/slice/Algorithm/AlgorithmAction'
import { isAlgoChild, isAlgoParent } from 'store/slice/Algorithm/AlgorithmUtils'
import { NODE_DATA_TYPE, NODE_DATA_TYPE_SET } from 'const/NodeData'

const useStyles = makeStyles({
  root: {
    // height: 240,
    flexGrow: 1,
    background: '#F0F0F0',
    height: 1000,
  },
})

export const SideBar = React.memo(() => {
  const dispatch = useDispatch()
  const classes = useStyles()
  const algoList = useSelector(algoListSelector, algoListEqualityFn)

  useEffect(() => {
    if (Object.keys(algoList).length === 0) {
      dispatch(getAlgoList())
    }
  }, [dispatch, algoList])

  const onDragStart = (
    event: DragEvent,
    nodeName: string,
    nodeDataType: NODE_DATA_TYPE,
    path?: string,
  ) => {
    if (event.dataTransfer != null) {
      event.dataTransfer.setData('nodeName', nodeName)
      event.dataTransfer.setData('type', nodeDataType)
      if (path != null) {
        event.dataTransfer.setData('path', path)
      }
      event.dataTransfer.effectAllowed = 'move'
    }
  }

  return (
    <TreeView
      className={classes.root}
      defaultCollapseIcon={<ExpandMoreIcon />}
      defaultExpandIcon={<ChevronRightIcon />}
    >
      <TreeItem nodeId="Data" label="Data">
        <TreeItem
          nodeId="image"
          label="image"
          onDragStart={(event: DragEvent) =>
            onDragStart(event, 'image', NODE_DATA_TYPE_SET.IMAGE)
          }
          draggable
        />
        <TreeItem
          nodeId="csv"
          label="csv"
          onDragStart={(event: DragEvent) =>
            onDragStart(event, 'csv', NODE_DATA_TYPE_SET.CSV)
          }
          draggable
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponent
            name={name}
            node={node}
            onDragStart={(event, nodeName, path) =>
              onDragStart(event, nodeName, NODE_DATA_TYPE_SET.ALGO, path)
            }
            key={i.toFixed()}
          />
        ))}
      </TreeItem>
    </TreeView>
  )
})

const AlgoNodeComponent = React.memo<{
  name: string
  node: AlgoNodeType
  onDragStart: (event: DragEvent, nodeName: string, path?: string) => void
}>(({ name, node, onDragStart }) => {
  if (isAlgoChild(node)) {
    return (
      <TreeItem
        nodeId={name}
        label={name}
        onDragStart={(event: DragEvent) => onDragStart(event, name, node.path)}
        draggable
      />
    )
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        {Object.entries(node.children).map(([name, node], i) => (
          <AlgoNodeComponent
            name={name}
            node={node}
            onDragStart={onDragStart}
            key={i.toFixed()}
          />
        ))}
      </TreeItem>
    )
  }
})

function algoListEqualityFn(a: AlgoListType, b: AlgoListType) {
  const aArray = Object.entries(a)
  const bArray = Object.entries(b)
  return (
    a === b ||
    (aArray.length === bArray.length &&
      aArray.every(([aKey, aValue], i) => {
        const [bKey, bValue] = bArray[i]
        return bKey === aKey && algoNodeEqualityFn(bValue, aValue)
      }))
  )
}

function algoNodeEqualityFn(a: AlgoNodeType, b: AlgoNodeType): boolean {
  if (a === b) {
    return true
  }
  if (isAlgoChild(a) && isAlgoChild(b)) {
    return arrayEqualityFn(a.args, b.args)
  } else if (isAlgoParent(a) && isAlgoParent(b)) {
    const aArray = Object.entries(a.children)
    const bArray = Object.entries(b.children)
    return (
      a.children === b.children ||
      (aArray.length === bArray.length &&
        aArray.every(([aKey, aValue], i) => {
          const [bKey, bValue] = bArray[i]
          return bKey === aKey && algoNodeEqualityFn(bValue, aValue)
        }))
    )
  } else {
    return false
  }
}

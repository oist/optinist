import React, { useEffect, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

import {
  selectAlgorithmListIsLated,
  selectAlgorithmListTree,
} from 'store/slice/AlgorithmList/AlgorithmListSelectors'
import { AlgorithmNodeType } from 'store/slice/AlgorithmList/AlgorithmListType'
import { isAlgoChild } from 'store/slice/AlgorithmList/AlgorithmListUtils'
import { getAlgoList } from 'store/slice/AlgorithmList/AlgorithmListActions'
import { NODE_TYPE_SET } from 'store/slice/FlowElement/FlowElementType'
import { FILE_TYPE, FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'

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
            onDataNodeDragStart(event, 'ImageData', FILE_TYPE_SET.IMAGE)
          }
          draggable
        />
        <TreeItem
          nodeId="csv"
          label="csv"
          onDragStart={(event: DragEvent) =>
            onDataNodeDragStart(event, 'CsvData', FILE_TYPE_SET.CSV)
          }
          draggable
        />
      </TreeItem>
      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.entries(algoList).map(([name, node], i) => (
          <AlgoNodeComponent
            name={name}
            node={node}
            onDragStart={(event, nodeName, functionPath) =>
              onAlgoNodeDragStart(event, nodeName, functionPath)
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
  node: AlgorithmNodeType
  onDragStart: (
    event: DragEvent,
    nodeName: string,
    functionPath: string,
  ) => void
}>(({ name, node, onDragStart }) => {
  if (isAlgoChild(node)) {
    return (
      <TreeItem
        nodeId={name}
        label={name}
        onDragStart={(event: DragEvent) =>
          onDragStart(event, name, node.functionPath)
        }
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

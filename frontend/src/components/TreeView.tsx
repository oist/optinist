import React, { useEffect, DragEvent } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

import { algoListSelector } from 'redux/slice/Algorithm/AlgorithmSelector'
import { AlgoListType } from 'redux/slice/Algorithm/AlgorithmType'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { getAlgoList } from 'redux/slice/Algorithm/AlgorithmAction'

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

  const onDragStart = (event: DragEvent, nodeName: string) => {
    if (event.dataTransfer != null) {
      event.dataTransfer.setData('application/reactflow', nodeName)
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
          nodeId="data"
          label="data"
          onDragStart={(event: DragEvent) => onDragStart(event, 'data')}
          draggable
        />
      </TreeItem>

      <TreeItem nodeId="Algorithm" label="Algorithm">
        {Object.keys(algoList).map((name) => (
          <TreeItem
            key={name}
            nodeId={name}
            label={name}
            onDragStart={(event: DragEvent) => onDragStart(event, name)}
            draggable
          />
        ))}
      </TreeItem>
    </TreeView>
  )
})

function algoListEqualityFn(a: AlgoListType, b: AlgoListType) {
  const aArray = Object.entries(a)
  const bArray = Object.entries(b)
  return (
    a === b ||
    (aArray.length === bArray.length &&
      aArray.every(([aKey, aValue], i) => {
        const [bKey, bValue] = bArray[i]
        return bKey === aKey && arrayEqualityFn(bValue.args, aValue.args)
      }))
  )
}

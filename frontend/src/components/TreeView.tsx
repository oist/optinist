import React, { useState, useEffect, DragEvent } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

import axios from 'axios'

const useStyles = makeStyles({
  root: {
    // height: 240,
    flexGrow: 1,
    background: '#F0F0F0',
    height: 1000,
  },
})

export const SideBar = React.memo(() => {
  const classes = useStyles()
  const [algoList, setAlgoList] = useState([])

  useEffect(() => {
    axios.get('http://localhost:8000/algolist').then((res) => {
      setAlgoList(res.data)
    })
  }, [])

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
        {algoList.map((name) => (
          <TreeItem
            key={name}
            nodeId={name}
            label={name}
            onDragStart={(event: DragEvent) => onDragStart(event, name)}
            draggable
          />
        ))}
      </TreeItem>

      <TreeItem nodeId="Output" label="Output">
        <TreeItem
          nodeId="output"
          label="output"
          onDragStart={(event: DragEvent) => onDragStart(event, 'output')}
          draggable
        />
      </TreeItem>
    </TreeView>
  )
})

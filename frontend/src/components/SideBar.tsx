import { DragEvent } from 'react'
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

const useStyles = makeStyles({
  root: {
    // height: 240,
    flexGrow: 1,
    background: '#F0F0F0',
    height: 1000,
  },
})

const SideBar = () => {
  const classes = useStyles()

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
          nodeId="data1"
          label="data1"
          onDragStart={(event: DragEvent) => onDragStart(event, 'data1')}
          draggable
        />
        <TreeItem
          nodeId="data2"
          label="data2"
          onDragStart={(event: DragEvent) => onDragStart(event, 'data2')}
          draggable
        />
        <TreeItem
          nodeId="data3"
          label="data3"
          onDragStart={(event: DragEvent) => onDragStart(event, 'data3')}
          draggable
        />
      </TreeItem>

      <TreeItem nodeId="Algorithm" label="Algorithm">
        <TreeItem
          nodeId="algo1"
          label="algo1"
          onDragStart={(event: DragEvent) => onDragStart(event, 'algo1')}
          draggable
        />
        <TreeItem
          nodeId="algo2"
          label="algo2"
          onDragStart={(event: DragEvent) => onDragStart(event, 'algo2')}
          draggable
        />
        <TreeItem
          nodeId="algo3"
          label="algo3"
          onDragStart={(event: DragEvent) => onDragStart(event, 'algo3')}
          draggable
        />
      </TreeItem>

      <TreeItem nodeId="output" label="Output">
        <TreeItem
          nodeId="output1"
          label="output1"
          onDragStart={(event: DragEvent) => onDragStart(event, 'output1')}
          draggable
        />
        <TreeItem
          nodeId="output2"
          label="output2"
          onDragStart={(event: DragEvent) => onDragStart(event, 'output3')}
          draggable
        />
        <TreeItem
          nodeId="output3"
          label="output3"
          onDragStart={(event: DragEvent) => onDragStart(event, 'output3')}
          draggable
        />
      </TreeItem>
    </TreeView>
  )
}

export default SideBar

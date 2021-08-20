// import { DragEvent } from 'react';
import { makeStyles } from '@material-ui/core/styles'
import TreeView from '@material-ui/lab/TreeView'
import TreeItem from '@material-ui/lab/TreeItem'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

// const onDragStart = (event: DragEvent, nodeType: string) => {
//   event.dataTransfer.setData('application/reactflow', nodeType);
//   event.dataTransfer.effectAllowed = 'move';
// };

const useStyles = makeStyles({
  root: {
    height: 240,
    flexGrow: 1,
    maxWidth: 400,
  },
})

const Sidebar = () => {
  const classes = useStyles()

  return (
    <TreeView
      className={classes.root}
      defaultCollapseIcon={<ExpandMoreIcon />}
      defaultExpandIcon={<ChevronRightIcon />}
    >
      <TreeItem nodeId="1" label="Applications">
        <TreeItem nodeId="2" label="Calendar" />
        <TreeItem nodeId="3" label="Chrome" />
        <TreeItem nodeId="4" label="Webstorm" />
      </TreeItem>
    </TreeView>
    // <aside>
    //   <p
    //     onDragStart={(event: DragEvent) => onDragStart(event, 'input')} draggable>
    //     Input Node
    //   </p>
    //   <p
    //     onDragStart={(event: DragEvent) => onDragStart(event, 'default')} draggable>
    //     Default Node
    //   </p>
    // </aside>
  )
}

export default Sidebar

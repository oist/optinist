import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Popover from '@material-ui/core/Popover'
import Typography from '@material-ui/core/Typography'
import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import Button from '@material-ui/core/Button'
import { nwbListSelector } from 'store/slice/NWB/NWBSelector'
import { NWBNodeType } from 'store/slice/NWB/NWBType'
import { getNWBParams } from 'store/slice/NWB/NWBAction'
import 'style/nwb.css'
import { TreeView, TreeItem } from '@material-ui/lab'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'

export const NWB = React.memo(() => {
  const [anchorEl, setAnchorEl] = React.useState<HTMLButtonElement | null>(null)

  const handleClick = (event: React.MouseEvent<HTMLButtonElement>) => {
    setAnchorEl(event.currentTarget)
  }

  const handleClose = () => {
    setAnchorEl(null)
  }

  const open = Boolean(anchorEl)
  const id = open ? 'simple-popover' : undefined

  return (
    <div>
      <Button
        aria-describedby={id}
        variant="contained"
        color="default"
        onClick={handleClick}
      >
        NWB setting
      </Button>
      <Popover
        id={id}
        open={open}
        anchorEl={anchorEl}
        onClose={handleClose}
        anchorOrigin={{
          vertical: 'bottom',
          horizontal: 'left',
        }}
        transformOrigin={{
          vertical: 'top',
          horizontal: 'center',
        }}
      >
        <ParamBar />
      </Popover>
    </div>
  )
})

const ParamBar = React.memo(() => {
  const nwbList = useSelector(nwbListSelector)

  const dispatch = useDispatch()

  useEffect(() => {
    if (Object.keys(nwbList).length === 0) {
      dispatch(getNWBParams())
    }
  })
  return (
    <div className="nwbParam">
      <h3>NWB Parameter</h3>
      <TreeView
        defaultCollapseIcon={<ExpandMoreIcon />}
        defaultExpandIcon={<ChevronRightIcon />}
      >
        {Object.entries(nwbList).map(([name, node], i) => (
          <ParamComponent key={name} name={name} node={node} />
        ))}
      </TreeView>
    </div>
  )
})

const ParamComponent = React.memo<{
  name: string
  node: NWBNodeType
}>(({ name, node }) => {
  if (node.type == 'child') {
    return <ParamItem name={name} value={node.value} />
  } else {
    return (
      <TreeItem nodeId={name} label={name}>
        <div key={name}>
          {Object.entries(node.children).map(([name, childNode], i) => (
            <ParamComponent name={name} node={childNode} />
          ))}
        </div>
      </TreeItem>
    )
  }
})

const ParamItem = React.memo<{
  name: string
  value: unknown
}>(({ name, value }) => {
  return (
    <Box
      sx={{
        display: 'flex',
        marginTop: 16,
        marginBottom: 16,
        alignItems: 'center',
      }}
    >
      <Box
        style={{ verticalAlign: 'middle' }}
        sx={{
          flexGrow: 1,
          width: '50%',
        }}
      >
        <Typography>{name}</Typography>
      </Box>
      <Box sx={{ width: '50%' }}>
        <TextField
          value={String(value)}
          // onChange={onChange}
        />
      </Box>
    </Box>
  )
})

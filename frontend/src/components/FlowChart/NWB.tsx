import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'

import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import Button from '@material-ui/core/Button'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import MuiAccordion from '@material-ui/core/Accordion'
import AccordionDetails from '@material-ui/core/AccordionDetails'
import AccordionSummary from '@material-ui/core/AccordionSummary'
import Typography from '@material-ui/core/Typography'
import { withStyles } from '@material-ui/core/styles'
import { selectNwbList } from 'store/slice/NWB/NWBSelectors'
import { NWBChild, NWBNodeType } from 'store/slice/NWB/NWBType'
import { updateParam } from 'store/slice/NWB/NWBSlice'
import { getNWBParams } from 'store/slice/NWB/NWBAction'
import { toggleNwb } from 'store/slice/RightDrawer/RightDrawerSlice'

export const NWBSettingButton = React.memo(() => {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleNwb())
  }
  return (
    <Button variant="contained" color="default" onClick={handleClick}>
      NWB setting
    </Button>
  )
})

export const NWBSettingContents = React.memo(() => {
  const nwbList = useSelector(selectNwbList)

  const dispatch = useDispatch()

  useEffect(() => {
    if (Object.keys(nwbList).length === 0) {
      dispatch(getNWBParams())
    }
  })
  return (
    <div className="nwbParam" style={{ margin: 4 }}>
      {Object.entries(nwbList).map(([name, node], i) => (
        <ParamComponent key={name} name={name} node={node} />
      ))}
    </div>
  )
})

const ParamComponent = React.memo<{
  name: string
  node: NWBNodeType
}>(({ name, node }) => {
  if (node.type === 'child') {
    return <ParamItem name={name} node={node} />
  } else {
    return (
      <Accordion square>
        <AccordionSummary
          expandIcon={<ExpandMoreIcon />}
          aria-controls="panel1bh-content"
          id="panel1bh-header"
        >
          {name}
        </AccordionSummary>
        <AccordionDetails>
          <div>
            {Object.entries(node.children).map(([name, childNode], i) => (
              <ParamComponent name={name} node={childNode} />
            ))}
          </div>
        </AccordionDetails>
      </Accordion>
    )
  }
})

const ParamItem = React.memo<{
  name: string
  node: NWBChild
}>(({ name, node: { value, path } }) => {
  const dispatch = useDispatch()
  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    dispatch(updateParam({ paramPath: path, newValue: event.target.value }))
  }
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
        <Typography style={{ overflow: 'scroll' }}>{name}</Typography>
      </Box>
      <Box sx={{ width: '50%' }}>
        <TextField value={String(value)} onChange={onChange} />
      </Box>
    </Box>
  )
})

const Accordion = withStyles((theme) => ({
  root: {
    border: '1px solid',
    borderColor: theme.palette.divider,
    boxShadow: 'none',
    '&:not(:last-child)': {
      borderBottom: 0,
    },
    '&:before': {
      display: 'none',
    },
    '&$expanded': {
      margin: 'auto',
    },
  },
  expanded: {},
}))(MuiAccordion)

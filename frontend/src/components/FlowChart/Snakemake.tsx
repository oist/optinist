import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'

import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import Button from '@material-ui/core/Button'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import AccordionDetails from '@material-ui/core/AccordionDetails'
import AccordionSummary from '@material-ui/core/AccordionSummary'
import Typography from '@material-ui/core/Typography'

import { selectSnakemakeList } from 'store/slice/Snakemake/SnakemakeSelectors'
import {
  SnakemakeChild,
  SnakemakeNodeType,
} from 'store/slice/Snakemake/SnakemakeType'
import { updateParam } from 'store/slice/Snakemake/SnakemakeSlice'
import { getSnakemakeParams } from 'store/slice/Snakemake/SnakemakeAction'
import { toggleSnakemake } from 'store/slice/RightDrawer/RightDrawerSlice'
import { Accordion } from 'components/Accordion'

export const SnakemakeButton = React.memo(() => {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleSnakemake())
  }
  return (
    <Button variant="contained" color="default" onClick={handleClick}>
      Snakemake
    </Button>
  )
})

export const SnakemakeContents = React.memo(() => {
  const SnakemakeList = useSelector(selectSnakemakeList)
  console.log(SnakemakeList)

  const dispatch = useDispatch()

  useEffect(() => {
    if (Object.keys(SnakemakeList).length === 0) {
      dispatch(getSnakemakeParams())
    }
  })
  return (
    <div className="SnakemakeParam" style={{ margin: 4 }}>
      {Object.entries(SnakemakeList).map(([name, node], i) => (
        <ParamComponent key={name} name={name} node={node} />
      ))}
    </div>
  )
})

const ParamComponent = React.memo<{
  name: string
  node: SnakemakeNodeType
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
  node: SnakemakeChild
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

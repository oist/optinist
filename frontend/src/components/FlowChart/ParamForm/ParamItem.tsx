import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { AnyAction } from '@reduxjs/toolkit'
import Switch from '@material-ui/core/Switch'
import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import ExpandMoreIcon from '@material-ui/icons/ExpandMore'
import AccordionDetails from '@material-ui/core/AccordionDetails'
import AccordionSummary from '@material-ui/core/AccordionSummary'
import Typography from '@material-ui/core/Typography'

import { updateParam } from 'store/slice/AlgorithmNode/AlgorithmNodeSlice'
import {
  selectAlgorithmParamsValue,
  selectAlgorithmParam,
} from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { isParamChild } from 'store/slice/AlgorithmNode/AlgorithmNodeUtils'
import { ParamType } from 'store/slice/AlgorithmNode/AlgorithmNodeType'
import { ParamFormContext } from '../RightDrawer'
import { Accordion } from 'components/Accordion'

export type ParamItemContainerProps = {
  paramKey: string
}

export const ParamItemContainer = React.memo<{ paramKey: string }>(
  ({ paramKey }) => {
    const nodeId = React.useContext(ParamFormContext)
    const param = useSelector(selectAlgorithmParam(nodeId, paramKey)) // 一階層目
    return <ParamItem paramKey={paramKey} param={param} />
  },
)

type ParamItemProps = {
  paramKey: string
  param: ParamType
}

const ParamItem = React.memo<ParamItemProps>(({ paramKey, param }) => {
  if (isParamChild(param)) {
    return <ParamChildItem path={param.path} name={paramKey} />
  } else {
    return (
      <Accordion square>
        <AccordionSummary expandIcon={<ExpandMoreIcon />}>
          {paramKey}
        </AccordionSummary>
        <AccordionDetails>
          <div>
            {Object.entries(param.children).map(([paramKey, param], i) => (
              <ParamItem param={param} paramKey={paramKey} />
            ))}
          </div>
        </AccordionDetails>
      </Accordion>
    )
  }
})

type ParamChildItemProps = {
  path: string
}

const ParamChildItem = React.memo<ParamChildItemProps & { name: string }>(
  ({ path, name }) => {
    return (
      <Box
        sx={{
          display: 'flex',
          marginTop: 16,
          marginBottom: 16,
          alignItems: 'center',
          overflow: 'scroll',
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
          <ParamItemForValueType path={path} />
        </Box>
      </Box>
    )
  },
)

const ParamItemForValueType = React.memo<ParamChildItemProps>(({ path }) => {
  const [value] = useParamValueUpdate(path)
  if (typeof value === 'number') {
    return <ParamItemForNumber path={path} />
  } else if (typeof value === 'string') {
    return <ParamItemForString path={path} />
  } else if (typeof value === 'boolean') {
    return <ParamItemForBoolean path={path} />
  } else {
    return <ParamItemForString path={path} />
  }
})

const ParamItemForString = React.memo<ParamChildItemProps>(({ path }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(path) // string値以外も想定されるため
  const onChange = (
    e: React.ChangeEvent<HTMLTextAreaElement | HTMLInputElement>,
  ) => {
    const newValue = e.target.value as string
    dispatch(updateParamAction(newValue))
  }
  return <TextField value={value} onChange={onChange} multiline />
})

const ParamItemForNumber = React.memo<ParamChildItemProps>(({ path }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(path)
  if (typeof value === 'number') {
    const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
      const newValue =
        event.target.value === '' ? '' : Number(event.target.value)
      if (typeof newValue === 'number') {
        dispatch(updateParamAction(newValue))
      }
    }
    return (
      <TextField
        type="number"
        InputLabelProps={{
          shrink: true,
        }}
        value={value}
        onChange={onChange}
      />
    )
  } else {
    return null
  }
})

const ParamItemForBoolean = React.memo<ParamChildItemProps>(({ path }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(path)
  if (typeof value === 'boolean') {
    const onChange = () => {
      dispatch(updateParamAction(!value))
    }
    return <Switch checked={value} onChange={onChange} />
  } else {
    return null
  }
})

function useParamValueUpdate(
  path: string,
): [unknown, (newValue: unknown) => AnyAction] {
  const nodeId = React.useContext(ParamFormContext)
  const value = useSelector(selectAlgorithmParamsValue(nodeId, path))
  const updateParamAction = (newValue: unknown) => {
    return updateParam({ nodeId, path, newValue })
  }
  return [value, updateParamAction]
}

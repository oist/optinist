import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { AnyAction } from '@reduxjs/toolkit'
import Switch from '@material-ui/core/Switch'
import TextField from '@material-ui/core/TextField'
import Box from '@material-ui/core/Box'
import Typography from '@material-ui/core/Typography'

import { updateParam } from 'store/slice/AlgorithmNode/AlgorithmNodeSlice'
import { selectAlgorithmParamsValue } from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { ParamFormTabContext } from 'App'

type ParamItemProps = {
  paramKey: string
}

export const ParamItemContainer = React.memo<ParamItemProps>(({ paramKey }) => {
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
        <Typography>{paramKey}</Typography>
      </Box>
      <Box sx={{ width: '50%' }}>
        <ParamItemForValueType paramKey={paramKey} />
      </Box>
    </Box>
  )
})

const ParamItemForValueType = React.memo<ParamItemProps>(({ paramKey }) => {
  const [value] = useParamValueUpdate(paramKey)
  if (typeof value === 'number') {
    return <ParamItemForNumber paramKey={paramKey} />
  } else if (typeof value === 'string') {
    return <ParamItemForString paramKey={paramKey} />
  } else if (typeof value === 'boolean') {
    return <ParamItemForBoolean paramKey={paramKey} />
  } else {
    return <ParamItemForString paramKey={paramKey} />
  }
})

const ParamItemForString = React.memo<ParamItemProps>(({ paramKey }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(paramKey) // string値以外も想定されるため
  const onChange = (
    e: React.ChangeEvent<HTMLTextAreaElement | HTMLInputElement>,
  ) => {
    const newValue = e.target.value as string
    dispatch(updateParamAction(newValue))
  }
  return <TextField value={value} onChange={onChange} multiline />
})

const ParamItemForNumber = React.memo<ParamItemProps>(({ paramKey }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(paramKey)
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

const ParamItemForBoolean = React.memo<ParamItemProps>(({ paramKey }) => {
  const dispatch = useDispatch()
  const [value, updateParamAction] = useParamValueUpdate(paramKey)
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
  paramKey: string,
): [unknown, (newValue: unknown) => AnyAction] {
  const nodeId = React.useContext(ParamFormTabContext)
  const value = useSelector(selectAlgorithmParamsValue(nodeId, paramKey))
  const updateParamAction = (newValue: unknown) => {
    return updateParam({ nodeId, paramKey, newValue })
  }
  return [value, updateParamAction]
}

import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { AnyAction } from '@reduxjs/toolkit'
import Switch from '@mui/material/Switch'
import TextField from '@mui/material/TextField'
import Box from '@mui/material/Box'
import ExpandMoreIcon from '@mui/icons-material/ExpandMore'
import AccordionDetails from '@mui/material/AccordionDetails'
import AccordionSummary from '@mui/material/AccordionSummary'
import Typography from '@mui/material/Typography'

import { isParamChild } from 'utils/param/ParamUtils'
import { ParamType } from 'utils/param/ParamType'
import { RootState } from 'store/store'
import { Accordion } from 'components/Accordion'

type ParamSelectorType = (paramKey: string) => (state: RootState) => ParamType
type ParamValueSelectorType = (path: string) => (state: RootState) => unknown
type ParamUpdateActionCreatorType = (
  path: string,
  newValue: unknown,
) => AnyAction

export type CreateParamFormItemComponentProps = {
  paramSelector: ParamSelectorType
  paramValueSelector: ParamValueSelectorType
  paramUpdateActionCreator: ParamUpdateActionCreatorType
}

export function createParamFormItemComponent({
  paramSelector,
  paramValueSelector,
  paramUpdateActionCreator,
}: CreateParamFormItemComponentProps): React.FC<{ paramKey: string }> {
  function useParamValueUpdate(
    path: string,
  ): [unknown, (newValue: unknown) => AnyAction] {
    const value = useSelector(paramValueSelector(path))
    const updateParamAction = (newValue: unknown) => {
      return paramUpdateActionCreator(path, newValue)
    }
    return [value, updateParamAction]
  }
  const ParamItemForString = React.memo<ParamChildItemProps>(({ path }) => {
    const dispatch = useDispatch()
    const [value, updateParamAction] = useParamValueUpdate(path)
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
  const ParamChildItem = React.memo<ParamChildItemProps & { name: string }>(
    ({ path, name }) => {
      return (
        <Box
          sx={{
            display: 'flex',
            marginTop: (theme) => theme.spacing(2),
            marginBottom: (theme) => theme.spacing(2),
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
  const ParamItem = React.memo<ParamItemProps>(({ paramKey, param }) => {
    if (isParamChild(param)) {
      return <ParamChildItem path={param.path} name={paramKey} />
    } else {
      return (
        <Accordion>
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
  return React.memo<{ paramKey: string }>(({ paramKey }) => {
    const param = useSelector(paramSelector(paramKey)) // 一階層目
    return <ParamItem paramKey={paramKey} param={param} />
  })
}

type ParamItemProps = {
  paramKey: string
  param: ParamType
}

type ParamChildItemProps = {
  path: string
}

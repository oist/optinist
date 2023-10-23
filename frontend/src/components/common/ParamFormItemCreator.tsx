import { ChangeEvent, FC, memo, useRef } from "react"
import { useDispatch, useSelector } from "react-redux"

import ExpandMoreIcon from "@mui/icons-material/ExpandMore"
import AccordionDetails from "@mui/material/AccordionDetails"
import AccordionSummary from "@mui/material/AccordionSummary"
import Box from "@mui/material/Box"
import Switch from "@mui/material/Switch"
import TextField from "@mui/material/TextField"
import Typography from "@mui/material/Typography"
import { AnyAction } from "@reduxjs/toolkit"

import { Accordion } from "components/common/Accordion"
import { RootState } from "store/store"
import { ParamType } from "utils/param/ParamType"
import { isParamChild } from "utils/param/ParamUtils"

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
}: CreateParamFormItemComponentProps): FC<{ paramKey: string }> {
  function useParamValueUpdate(
    path: string,
  ): [unknown, (newValue: unknown) => AnyAction] {
    const value = useSelector(paramValueSelector(path))
    const updateParamAction = (newValue: unknown) => {
      return paramUpdateActionCreator(path, newValue)
    }
    return [value, updateParamAction]
  }
  const ParamItemForString = memo(function ParamItemForString({
    path,
  }: ParamChildItemProps) {
    const dispatch = useDispatch()
    const [value, updateParamAction] = useParamValueUpdate(path)
    const isArray = useRef(Array.isArray(value))
    const onChange = (
      e: ChangeEvent<HTMLTextAreaElement | HTMLInputElement>,
    ) => {
      const newValue = e.target.value as string
      dispatch(updateParamAction(newValue))
    }

    const onBlur = (e: ChangeEvent<HTMLTextAreaElement | HTMLInputElement>) => {
      const newValue = e.target.value as string
      dispatch(
        updateParamAction(
          newValue
            .split(",")
            .filter(Boolean)
            .map((e) => Number(e)),
        ),
      )
    }
    return (
      <TextField
        value={value}
        onChange={onChange}
        multiline
        onBlur={isArray ? onBlur : undefined}
      />
    )
  })
  const ParamItemForNumber = memo(function ParamItemForNumber({
    path,
  }: ParamChildItemProps) {
    const dispatch = useDispatch()
    const [value, updateParamAction] = useParamValueUpdate(path)
    if (typeof value === "number") {
      const onChange = (event: ChangeEvent<HTMLInputElement>) => {
        const newValue =
          event.target.value === "" ? "" : Number(event.target.value)
        if (typeof newValue === "number") {
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
  const ParamItemForBoolean = memo(function ParamItemForBoolean({
    path,
  }: ParamChildItemProps) {
    const dispatch = useDispatch()
    const [value, updateParamAction] = useParamValueUpdate(path)
    if (typeof value === "boolean") {
      const onChange = () => {
        dispatch(updateParamAction(!value))
      }
      return <Switch checked={value} onChange={onChange} />
    } else {
      return null
    }
  })
  const ParamItemForValueType = memo(function ParamItemForValueType({
    path,
  }: ParamChildItemProps) {
    const [value] = useParamValueUpdate(path)
    if (typeof value === "number") {
      return <ParamItemForNumber path={path} />
    } else if (typeof value === "string") {
      return <ParamItemForString path={path} />
    } else if (typeof value === "boolean") {
      return <ParamItemForBoolean path={path} />
    } else {
      return <ParamItemForString path={path} />
    }
  })
  const ParamChildItem = memo(function ParamChildItem({
    path,
    name,
  }: ParamChildItemWithNameProps) {
    return (
      <Box
        sx={{
          display: "flex",
          marginTop: (theme) => theme.spacing(2),
          marginBottom: (theme) => theme.spacing(2),
          alignItems: "center",
          overflow: "scroll",
        }}
      >
        <Box
          style={{ verticalAlign: "middle" }}
          sx={{
            flexGrow: 1,
            width: "50%",
          }}
        >
          <Typography style={{ overflow: "scroll" }}>{name}</Typography>
        </Box>
        <Box sx={{ width: "50%" }}>
          <ParamItemForValueType path={path} />
        </Box>
      </Box>
    )
  })
  const ParamItem = memo(function ParamItem({
    paramKey,
    param,
  }: ParamItemProps) {
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
              {Object.entries(param.children).map(([paramKey, param]) => (
                <ParamItem key={paramKey} param={param} paramKey={paramKey} />
              ))}
            </div>
          </AccordionDetails>
        </Accordion>
      )
    }
  })
  return memo(function CreateParamFormItemComponen({
    paramKey,
  }: ParamKeyProps) {
    const param = useSelector(paramSelector(paramKey)) // 一階層目
    return <ParamItem paramKey={paramKey} param={param} />
  })
}

interface ParamKeyProps {
  paramKey: string
}

interface ParamItemProps extends ParamKeyProps {
  param: ParamType
}

interface ParamChildItemProps {
  path: string
}

interface ParamChildItemWithNameProps extends ParamChildItemProps {
  name: string
}

import { ChangeEvent, FC, memo, useContext, useRef } from "react"
import { useDispatch, useSelector } from "react-redux"

import ExpandMoreIcon from "@mui/icons-material/ExpandMore"
import { Tooltip } from "@mui/material"
import AccordionDetails from "@mui/material/AccordionDetails"
import AccordionSummary from "@mui/material/AccordionSummary"
import Box from "@mui/material/Box"
import Switch from "@mui/material/Switch"
import TextField from "@mui/material/TextField"
import Typography from "@mui/material/Typography"
import { AnyAction } from "@reduxjs/toolkit"

import { Accordion } from "components/common/Accordion"
import { DialogContext } from "components/Workspace/FlowChart/Dialog/DialogContext"
import { selectPipelineLatestUid } from "store/slice/Pipeline/PipelineSelectors"
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
  requireConfirm?: boolean
}

export function createParamFormItemComponent({
  paramSelector,
  paramValueSelector,
  paramUpdateActionCreator,
  requireConfirm,
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
    const currentWorkflowId = useSelector(selectPipelineLatestUid)
    const { onOpenClearWorkflowIdDialog } = useContext(DialogContext)

    const onChange = (
      e: ChangeEvent<HTMLTextAreaElement | HTMLInputElement>,
    ) => {
      const newValue = e.target.value as string
      if (requireConfirm && currentWorkflowId != null) {
        onOpenClearWorkflowIdDialog({
          open: true,
          handleOk: () => {
            dispatch(updateParamAction(newValue))
          },
          handleCancel: () => null,
        })
      } else {
        dispatch(updateParamAction(newValue))
      }
    }

    const splitValue = (value: string) =>
      value
        .split(",")
        .filter(Boolean)
        .map((e) => Number(e))

    const onBlur = (e: ChangeEvent<HTMLTextAreaElement | HTMLInputElement>) => {
      const newValue = e.target.value as string
      if (requireConfirm && currentWorkflowId != null) {
        onOpenClearWorkflowIdDialog({
          open: true,
          handleOk: () => {
            dispatch(updateParamAction(splitValue(newValue)))
          },
          handleCancel: () => null,
        })
      } else {
        dispatch(updateParamAction(splitValue(newValue)))
      }
    }
    const valueField = Array.isArray(value)
      ? value.toLocaleString()
      : (value as string) || ""
    return (
      <TextField
        value={valueField === undefined ? "" : valueField}
        onChange={onChange}
        multiline
        onBlur={isArray.current ? onBlur : undefined}
      />
    )
  })

  const ParamItemForNumber = memo(function ParamItemForNumber({
    path,
  }: ParamChildItemProps) {
    const dispatch = useDispatch()
    const [value, updateParamAction] = useParamValueUpdate(path)
    const currentWorkflowId = useSelector(selectPipelineLatestUid)
    const { onOpenClearWorkflowIdDialog } = useContext(DialogContext)

    if (typeof value === "number") {
      const onChange = (event: ChangeEvent<HTMLInputElement>) => {
        const newValue =
          event.target.value === "" ? "" : Number(event.target.value)
        if (typeof newValue === "number") {
          if (requireConfirm && currentWorkflowId != null) {
            onOpenClearWorkflowIdDialog({
              open: true,
              handleOk: () => {
                dispatch(updateParamAction(newValue))
              },
              handleCancel: () => null,
            })
          } else {
            dispatch(updateParamAction(newValue))
          }
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
    const currentWorkflowId = useSelector(selectPipelineLatestUid)
    const { onOpenClearWorkflowIdDialog } = useContext(DialogContext)

    if (typeof value === "boolean") {
      const onChange = () => {
        if (requireConfirm && currentWorkflowId != null) {
          onOpenClearWorkflowIdDialog({
            open: true,
            handleOk: () => {
              dispatch(updateParamAction(!value))
            },
            handleCancel: () => null,
          })
        } else {
          dispatch(updateParamAction(!value))
        }
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
            cursor: "default",
          }}
        >
          <Tooltip
            title={<span style={{ fontSize: 16 }}>{name}</span>}
            placement={"top"}
          >
            <Typography
              style={{
                overflow: "hidden",
                width: "90%",
                textOverflow: "ellipsis",
              }}
            >
              {name}
            </Typography>
          </Tooltip>
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

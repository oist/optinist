import React, { useEffect } from "react"
import { useSelector, useDispatch } from "react-redux"

import RouteIcon from "@mui/icons-material/Route"
import { IconButton, Tooltip } from "@mui/material"

import { createParamFormItemComponent } from "components/common/ParamFormItemCreator"
import { toggleSnakemake } from "store/slice/RightDrawer/RightDrawerSlice"
import { getSnakemakeParams } from "store/slice/Snakemake/SnakemakeAction"
import {
  selectSnakemakeParam,
  selectSnakemakeParamsKeyList,
  selectSnakemakeParamsValue,
} from "store/slice/Snakemake/SnakemakeSelectors"
import { updateParam } from "store/slice/Snakemake/SnakemakeSlice"
import { AppDispatch } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"



export const SnakemakeButton = React.memo(() => {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleSnakemake())
  }
  return (
    <Tooltip title="Snakemake settings">
      <IconButton onClick={handleClick}>
        <RouteIcon color="primary" />
      </IconButton>
    </Tooltip>
  )
})

export const SnakemakeContents = React.memo(() => {
  const dispatch = useDispatch<AppDispatch>()
  const paramKeyList = useSelector(
    selectSnakemakeParamsKeyList,
    arrayEqualityFn,
  )
  useEffect(() => {
    if (paramKeyList.length === 0) {
      dispatch(getSnakemakeParams())
    }
  })
  return (
    <div className="SnakemakeParam" style={{ margin: 4 }}>
      {paramKeyList.map((paramKey, i) => (
        <ParamItem key={i} paramKey={paramKey} />
      ))}
    </div>
  )
})

const ParamItem = React.memo<{ paramKey: string }>(({ paramKey }) => {
  const Component = createParamFormItemComponent({
    paramSelector: selectSnakemakeParam,
    paramValueSelector: selectSnakemakeParamsValue,
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

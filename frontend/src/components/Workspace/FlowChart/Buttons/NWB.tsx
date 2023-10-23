import { memo, useEffect } from "react"
import { useSelector, useDispatch } from "react-redux"

import TuneIcon from "@mui/icons-material/Tune"
import { IconButton, Tooltip } from "@mui/material"

import { createParamFormItemComponent } from "components/common/ParamFormItemCreator"
import { getNWBParams } from "store/slice/NWB/NWBAction"
import {
  selectNwbParam,
  selectNwbParamsKeyList,
  selectNwbParamsValue,
} from "store/slice/NWB/NWBSelectors"
import { updateParam } from "store/slice/NWB/NWBSlice"
import { toggleNwb } from "store/slice/RightDrawer/RightDrawerSlice"
import { AppDispatch } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"

export const NWBSettingButton = memo(function NWBSettingButton() {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleNwb())
  }
  return (
    <Tooltip title="NWB settings">
      <IconButton onClick={handleClick}>
        <TuneIcon color="primary" />
      </IconButton>
    </Tooltip>
  )
})

export const NWBSettingContents = memo(function NWBSettingContents() {
  const dispatch = useDispatch<AppDispatch>()

  const paramKeyList = useSelector(selectNwbParamsKeyList, arrayEqualityFn)
  useEffect(() => {
    if (paramKeyList.length === 0) {
      dispatch(getNWBParams())
    }
  })

  return (
    <div className="nwbParam" style={{ margin: 4 }}>
      {paramKeyList.map((paramKey, i) => (
        <ParamItem key={i} paramKey={paramKey} />
      ))}
    </div>
  )
})

interface ParamItemProps {
  paramKey: string
}

const ParamItem = memo(function ParamItem({ paramKey }: ParamItemProps) {
  const Component = createParamFormItemComponent({
    paramSelector: selectNwbParam,
    paramValueSelector: selectNwbParamsValue,
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

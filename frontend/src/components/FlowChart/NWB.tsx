import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'

import Button from '@mui/material/Button'
import TuneIcon from '@mui/icons-material/Tune'
import { updateParam } from 'store/slice/NWB/NWBSlice'
import { getNWBParams } from 'store/slice/NWB/NWBAction'
import { toggleNwb } from 'store/slice/RightDrawer/RightDrawerSlice'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectNwbParam,
  selectNwbParamsKeyList,
  selectNwbParamsValue,
} from 'store/slice/NWB/NWBSelectors'
import { createParamFormItemComponent } from 'components/common/ParamFormItemCreator'
import { AppDispatch } from 'store/store'

export const NWBSettingButton = React.memo(() => {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleNwb())
  }
  return (
    <Button
      variant="outlined"
      onClick={handleClick}
      sx={{
        margin: (theme) => theme.spacing(1),
      }}
      endIcon={<TuneIcon />}
    >
      NWB setting
    </Button>
  )
})

export const NWBSettingContents = React.memo(() => {
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

const ParamItem = React.memo<{ paramKey: string }>(({ paramKey }) => {
  const Component = createParamFormItemComponent({
    paramSelector: selectNwbParam,
    paramValueSelector: selectNwbParamsValue,
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

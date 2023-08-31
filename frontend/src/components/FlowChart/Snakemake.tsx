import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Button from '@mui/material/Button'
import TuneIcon from '@mui/icons-material/Tune'

import { updateParam } from 'store/slice/Snakemake/SnakemakeSlice'
import { getSnakemakeParams } from 'store/slice/Snakemake/SnakemakeAction'
import { toggleSnakemake } from 'store/slice/RightDrawer/RightDrawerSlice'
import {
  selectSnakemakeParam,
  selectSnakemakeParamsKeyList,
  selectSnakemakeParamsValue,
} from 'store/slice/Snakemake/SnakemakeSelectors'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { createParamFormItemComponent } from 'components/common/ParamFormItemCreator'
import { AppDispatch } from 'store/store'

export const SnakemakeButton = React.memo(() => {
  const dispatch = useDispatch()
  const handleClick = () => {
    dispatch(toggleSnakemake())
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
      Snakemake
    </Button>
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

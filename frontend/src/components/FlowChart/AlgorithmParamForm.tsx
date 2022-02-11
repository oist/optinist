import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Typography from '@mui/material/Typography'

import {
  selectAlgorithmName,
  selectAlgorithmParamsExit,
  selectAlgorithmParamsKeyList,
  selectAlgorithmParamsValue,
  selectAlgorithmParam,
} from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { updateParam } from 'store/slice/AlgorithmNode/AlgorithmNodeSlice'
import { getAlgoParams } from 'store/slice/AlgorithmNode/AlgorithmNodeActions'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { createParamFormItemComponent } from '../ParamFormItemCreator'
import { ParamFormContext } from './RightDrawer'

export const AlgorithmParamForm = React.memo(() => {
  const nodeId = React.useContext(ParamFormContext)
  const dispatch = useDispatch()
  const algoName = useSelector(selectAlgorithmName(nodeId))
  const algoParamIsLoaded = useSelector(selectAlgorithmParamsExit(nodeId))
  const paramKeyList = useSelector(
    selectAlgorithmParamsKeyList(nodeId),
    arrayEqualityFn,
  )
  useEffect(() => {
    if (!algoParamIsLoaded) {
      dispatch(getAlgoParams({ nodeId, algoName }))
    }
  }, [dispatch, nodeId, algoName, algoParamIsLoaded])
  return (
    <div style={{ padding: 8 }}>
      <Typography variant="h5">{algoName}</Typography>
      <div style={{ paddingLeft: 8 }}>
        {paramKeyList.map((paramKey) => (
          <ParamItem key={paramKey} paramKey={paramKey} />
        ))}
      </div>
    </div>
  )
})

const ParamItem = React.memo<{ paramKey: string }>(({ paramKey }) => {
  const nodeId = React.useContext(ParamFormContext)
  const Component = createParamFormItemComponent({
    paramSelector: (paramKey) => selectAlgorithmParam(nodeId, paramKey),
    paramValueSelector: (path) => selectAlgorithmParamsValue(nodeId, path),
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ nodeId, path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

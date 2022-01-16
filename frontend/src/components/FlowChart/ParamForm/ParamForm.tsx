import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import Typography from '@material-ui/core/Typography'

import {
  selectAlgorithmName,
  selectAlgorithmParamsExit,
  selectAlgorithmParamsKeyList,
} from 'store/slice/AlgorithmNode/AlgorithmNodeSelectors'
import { getAlgoParams } from 'store/slice/AlgorithmNode/AlgorithmNodeActions'

import { arrayEqualityFn } from 'utils/EqualityUtils'
import { ParamItemContainer } from './ParamItem'
import { ParamFormContext } from '../RightDrawer'

const ParamForm = React.memo(() => {
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
      <Typography variant="h5">
        {algoName}({nodeId})
      </Typography>
      <div style={{ paddingLeft: 8 }}>
        {paramKeyList.map((paramKey) => (
          <ParamItemContainer key={paramKey} paramKey={paramKey} />
        ))}
      </div>
    </div>
  )
})

export default ParamForm

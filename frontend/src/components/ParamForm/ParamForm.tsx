import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { nodeByIdSelector } from 'redux/slice/Element/ElementSelector'

import { getAlgoParams } from 'redux/slice/Algorithm/AlgorithmAction'
import { ParamItemContainer } from './ParamItem'
import Typography from '@material-ui/core/Typography'
import { isAlgoNodeData } from 'utils/ElementUtils'
import {
  algoParamByIdSelector,
  currentAlgoIdSelector,
  currentAlgoNameSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'

export const ParamForm = React.memo(() => {
  const currentAlgoId = useSelector(currentAlgoIdSelector)
  const currentAlgoName = useSelector(currentAlgoNameSelector)
  const currentNode = useSelector(nodeByIdSelector(currentAlgoId))
  const algoParam = useSelector(algoParamByIdSelector(currentAlgoId))
  const dispatch = useDispatch()

  useEffect(() => {
    if (isAlgoNodeData(currentNode) && currentNode.data && !algoParam) {
      const algoName = currentNode.data.label
      dispatch(getAlgoParams({ id: currentAlgoId, algoName }))
    }
  }, [dispatch, currentAlgoId, algoParam, currentNode])

  if (algoParam === undefined) {
    return null
  }

  return (
    <div style={{ padding: 8 }}>
      <Typography variant="h5">
        {currentAlgoName}({currentAlgoId})
      </Typography>
      <div style={{ paddingLeft: 8 }}>
        {Object.keys(algoParam).map((paramName) => (
          <ParamItemContainer key={paramName} paramKey={paramName} />
        ))}
      </div>
    </div>
  )
})

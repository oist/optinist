import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { nodeByIdSelector } from 'redux/slice/Element/ElementSelector'

import { getAlgoParams } from 'redux/slice/Algorithm/AlgorithmAction'
import { ParamItemContainer } from './ParamItem'
import Typography from '@material-ui/core/Typography'
import { isAlgoNodeData } from 'utils/ElementUtils'
import {
  algoParamByIdSelector,
  algoNameByIdSelector,
} from 'redux/slice/Algorithm/AlgorithmSelector'
import { NodeIdContext } from 'App'

export const ParamForm = React.memo(() => {
  const nodeId = React.useContext(NodeIdContext)
  const currentAlgoName = useSelector(algoNameByIdSelector(nodeId))
  const currentNode = useSelector(nodeByIdSelector(nodeId))
  const algoParam = useSelector(algoParamByIdSelector(nodeId))
  const dispatch = useDispatch()

  useEffect(() => {
    if (isAlgoNodeData(currentNode) && currentNode.data && !algoParam) {
      const algoName = currentNode.data.label
      dispatch(getAlgoParams({ id: nodeId, algoName }))
    }
  }, [dispatch, nodeId, algoParam, currentNode])

  if (algoParam === undefined) {
    return null
  }

  return (
    <div style={{ padding: 8 }}>
      <Typography variant="h5">
        {currentAlgoName}({nodeId})
      </Typography>
      <div style={{ paddingLeft: 8 }}>
        {Object.keys(algoParam).map((paramName) => (
          <ParamItemContainer key={paramName} paramKey={paramName} />
        ))}
      </div>
    </div>
  )
})

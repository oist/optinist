import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import {
  currentAlgoIdSelector,
  elementByIdSelector,
  algoParamByIdSelector,
} from 'redux/slice/Element/ElementSelector'
import { getAlgoParams } from 'redux/slice/Element/ElementAction'
import { ParamItemContainer } from './ParamItem'
import Typography from '@material-ui/core/Typography'
import { isAlgoNodeData } from 'redux/slice/Element/ElementUtils'

export const ParamForm = React.memo(() => {
  const currentAlgoId = useSelector(currentAlgoIdSelector)
  const currentNode = useSelector(elementByIdSelector(currentAlgoId))
  const algoParam = useSelector(algoParamByIdSelector(currentAlgoId))
  const dispatch = useDispatch()

  useEffect(() => {
    if (isAlgoNodeData(currentNode) && currentNode.data && !algoParam) {
      const algoName = currentNode.data.label
      dispatch(getAlgoParams({ id: currentAlgoId, algoName }))
    }
  }, [currentAlgoId])

  if (algoParam === undefined) {
    return null
  }

  return (
    <div style={{ padding: 8 }}>
      <Typography variant="h5">
        {algoParam.name}({currentAlgoId})
      </Typography>
      <div style={{ paddingLeft: 8 }}>
        {Object.keys(algoParam.param).map((paramName) => (
          <ParamItemContainer key={paramName} paramKey={paramName} />
        ))}
      </div>
    </div>
  )
})

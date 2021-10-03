import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import {
  currentElementIdSelector,
  elementByIdSelector,
  algoParamByIdSelector,
} from 'redux/slice/Element/ElementSelector'
import { getAlgoParams } from 'redux/slice/Element/ElementAction'
import { ParamItemContainer } from './ParamItem'
import Typography from '@material-ui/core/Typography'

export const ParamForm = React.memo(() => {
  const currentElementId = useSelector(currentElementIdSelector)
  const currentElement = useSelector(elementByIdSelector(currentElementId))
  const algoParam = useSelector(algoParamByIdSelector(currentElementId))
  const dispatch = useDispatch()

  useEffect(() => {
    if (
      currentElement &&
      currentElement.data &&
      currentElement?.type === 'default' && // algoはElementのtypeが'default'、それ以外のelementはparamを持たないのでdispatchしない
      !algoParam
    ) {
      const algoName = currentElement.data.label
      dispatch(getAlgoParams({ id: currentElementId, algoName }))
    }
  }, [currentElementId])

  if (algoParam === undefined) {
    return null
  }

  return (
    <div style={{ padding: 8 }}>
      <Typography variant="h5">
        {algoParam.name}({currentElementId})
      </Typography>
      <div style={{ paddingLeft: 8 }}>
        {Object.keys(algoParam.param).map((paramName) => (
          <ParamItemContainer key={paramName} paramKey={paramName} />
        ))}
      </div>
    </div>
  )
})

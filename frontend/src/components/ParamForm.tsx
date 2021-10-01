import React, { useEffect } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import {
  currentElementSelector,
  algoParamsSelector,
} from 'redux/slice/Element/ElementSelector'
import ParamItem from './paramFormList/ParamItem'
import { getAlgoParams } from 'redux/slice/Element/ElementAction'

export const ParamForm = React.memo(() => {
  const currentElement = useSelector(currentElementSelector)
  const algoParams = useSelector(algoParamsSelector)
  const dispatch = useDispatch()
  console.log(currentElement)

  useEffect(() => {
    if (!Object.prototype.hasOwnProperty.call(algoParams, currentElement)) {
      dispatch(getAlgoParams(currentElement))
    }
  }, [currentElement])

  if (!Object.prototype.hasOwnProperty.call(algoParams, currentElement)) {
    return null
  }

  return (
    <div>
      <h2>{currentElement}</h2>
      <ul>
        {Object.keys(algoParams[currentElement]).map((name) => (
          <ParamItem key={name} name={name} />
        ))}
      </ul>
    </div>
  )
})

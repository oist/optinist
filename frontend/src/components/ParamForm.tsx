import { useSelector } from 'react-redux'
import {
  currentElementSelector,
  algoParamsSelector,
} from 'redux/slice/Element/ElementSelector'
import ParamItem from './ParamItem'

const ParamForm = () => {
  const currentElement = useSelector(currentElementSelector)
  const algoParams = useSelector(algoParamsSelector)

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
}

export default ParamForm

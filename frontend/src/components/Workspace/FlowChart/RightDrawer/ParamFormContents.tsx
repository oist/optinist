import { createContext, FC } from "react"
import { useSelector } from "react-redux"

import { AlgorithmParamForm } from "components/Workspace/FlowChart/RightDrawer/AlgorithmParamForm"
import { selectRightDrawerCurrentNodeId } from "store/slice/RightDrawer/RightDrawerSelectors"

export const ParamFormContext = createContext<string>("")

export const ParamFormContents: FC = () => {
  const nodeId = useSelector(selectRightDrawerCurrentNodeId)
  if (nodeId != null) {
    return (
      <ParamFormContext.Provider value={nodeId}>
        <AlgorithmParamForm />
      </ParamFormContext.Provider>
    )
  } else {
    return null
  }
}

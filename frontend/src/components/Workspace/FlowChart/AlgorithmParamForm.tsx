import { memo, useContext, useEffect } from "react"
import { useSelector, useDispatch } from "react-redux"

import Typography from "@mui/material/Typography"

import { createParamFormItemComponent } from "components/common/ParamFormItemCreator"
import { ParamFormContext } from "components/Workspace/FlowChart/RightDrawer"
import { getAlgoParams } from "store/slice/AlgorithmNode/AlgorithmNodeActions"
import {
  selectAlgorithmName,
  selectAlgorithmParamsExit,
  selectAlgorithmParamsKeyList,
  selectAlgorithmParamsValue,
  selectAlgorithmParam,
} from "store/slice/AlgorithmNode/AlgorithmNodeSelectors"
import { updateParam } from "store/slice/AlgorithmNode/AlgorithmNodeSlice"
import { ParamItemProps } from "store/slice/RightDrawer/RightDrawerType"
import { AppDispatch } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"

export const AlgorithmParamForm = memo(function AlgorithmParamForm() {
  const nodeId = useContext<string>(ParamFormContext)
  const dispatch = useDispatch<AppDispatch>()
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

const ParamItem = memo(function ParamItem({ paramKey }: ParamItemProps) {
  const nodeId = useContext(ParamFormContext)
  const Component = createParamFormItemComponent({
    paramSelector: (paramKey) => selectAlgorithmParam(nodeId, paramKey),
    paramValueSelector: (path) => selectAlgorithmParamsValue(nodeId, path),
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ nodeId, path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

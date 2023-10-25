import { memo, useEffect } from "react"
import { useSelector, useDispatch } from "react-redux"

import { createParamFormItemComponent } from "components/common/ParamFormItemCreator"
import { ParamItemProps } from "store/slice/RightDrawer/RightDrawerType"
import { getSnakemakeParams } from "store/slice/Snakemake/SnakemakeAction"
import {
  selectSnakemakeParam,
  selectSnakemakeParamsKeyList,
  selectSnakemakeParamsValue,
} from "store/slice/Snakemake/SnakemakeSelectors"
import { updateParam } from "store/slice/Snakemake/SnakemakeSlice"
import { AppDispatch } from "store/store"
import { arrayEqualityFn } from "utils/EqualityUtils"

export const SnakemakeSettingContents = memo(
  function SnakemakeSettingContents() {
    const dispatch = useDispatch<AppDispatch>()
    const paramKeyList = useSelector(
      selectSnakemakeParamsKeyList,
      arrayEqualityFn,
    )
    useEffect(() => {
      if (paramKeyList.length === 0) {
        dispatch(getSnakemakeParams())
      }
    })
    return (
      <div className="SnakemakeParam" style={{ margin: 4 }}>
        {paramKeyList.map((paramKey, i) => (
          <ParamItem key={i} paramKey={paramKey} />
        ))}
      </div>
    )
  },
)

const ParamItem = memo(function ParamItem({ paramKey }: ParamItemProps) {
  const Component = createParamFormItemComponent({
    paramSelector: selectSnakemakeParam,
    paramValueSelector: selectSnakemakeParamsValue,
    paramUpdateActionCreator: (path, newValue) =>
      updateParam({ path, newValue }),
  })
  return <Component paramKey={paramKey} />
})

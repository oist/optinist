import { getChildParam } from 'store/utils/param/ParamUtils'
import { RootState } from '../../store'

export const selectNwb = (state: RootState) => state.nwb

export const selectNwbParams = (state: RootState) => selectNwb(state).params

export const selectNwbParamsKeyList = (state: RootState) =>
  Object.keys(selectNwb(state).params)

export const selectNwbParam = (paramKey: string) => (state: RootState) =>
  selectNwb(state).params[paramKey]

export const selectNwbParamsValue = (path: string) => (state: RootState) => {
  const params = selectNwb(state).params
  if (params != null) {
    const target = getChildParam(path, params)
    return target?.value
  } else {
    throw new Error('Param is null')
  }
}

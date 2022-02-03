import { getChildParam } from 'store/utils/param/ParamUtils'
import { RootState } from '../../store'

const selectSnakemake = (state: RootState) => state.snakemake
export const selectSnakemakeParams = (state: RootState) =>
  selectSnakemake(state).params

export const selectSnakemakeParamsKeyList = (state: RootState) =>
  Object.keys(selectSnakemakeParams(state))

export const selectSnakemakeParam = (paramKey: string) => (state: RootState) =>
  selectSnakemakeParams(state)[paramKey]

export const selectSnakemakeParamsValue =
  (path: string) => (state: RootState) => {
    const params = selectSnakemakeParams(state)
    if (params != null) {
      const target = getChildParam(path, params)
      return target?.value
    } else {
      throw new Error('Param is null')
    }
  }

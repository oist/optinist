import { RootState } from '../../store'

export const selectSnakemake = (state: RootState) => state.snakemake

export const selectSnakemakeList = (state: RootState) => {
  return selectSnakemake(state).SnakemakeList
}

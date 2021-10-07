export const OUTPUT_SLICE_NAME = 'Output'

export type Output = {
  currentOutputId: string
  outputData: {
    [id: string]: OutputData[]
  }
}

export type OutputData = {
  x: number | string
  y: number
}

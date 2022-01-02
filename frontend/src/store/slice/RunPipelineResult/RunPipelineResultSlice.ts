import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { RunPipelineDTO } from 'api/Run/Run'

type RunPipeline = {
  result: RunPipelineDTO
}

const initialState: RunPipeline = {
  result: {
    status: 'ready',
    message: '',
    // outputPaths: {},
  },
  // result: {
  //   message: 'completed',
  //   status: 'completed',
  //   requestId: 'X2j_lBQd5M9qM8N-a5bDD',
  //   outputPaths: {
  //     '/tmp/optinist/copy_image1/copy_image1.tif': {
  //       images: {
  //         path: 'null',
  //         type: 'images',
  //         max_index: 3000,
  //       },
  //     },
  //     'dummy/dummy_image2image': {
  //       image2image: {
  //         path: '/tmp/optinist/dummy_image2image/image.json',
  //         type: 'images',
  //         max_index: 30,
  //       },
  //     },
  //     'dummy/dummy_image2time': {
  //       image2time: {
  //         path: '/tmp/optinist/dummy_image2time/timeseries.json',
  //         type: 'timeseries',
  //       },
  //     },
  //   },
  // },
}

export const runPipelineResultSlice = createSlice({
  name: 'runPipelineResult',
  initialState,
  reducers: {
    reflectResult: (state, action: PayloadAction<RunPipelineDTO>) => {
      state.result = action.payload
    },
  },
})

export const { reflectResult } = runPipelineResultSlice.actions

export default runPipelineResultSlice.reducer

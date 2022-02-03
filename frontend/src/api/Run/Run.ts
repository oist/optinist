import { createApi, fetchBaseQuery } from '@reduxjs/toolkit/query/react'
import { BASE_URL, WS_BASE_URL } from 'const/API'
import { Edge, Node } from 'react-flow-renderer'
import { ParamMap } from 'store/utils/param/ParamType'

export type RunPipelineDTO = {
  status: string
  message: string
  name?: string
  outputPaths?: OutputPathsDTO
  requestId?: string
}

export type OutputPathsDTO = {
  [path: string]: {
    [key: string]: {
      path: string
      type: string
      max_index?: number
    }
  }
}

export const webSocketApi = createApi({
  reducerPath: 'webSocketApi',
  baseQuery: fetchBaseQuery({ baseUrl: `${BASE_URL}/` }),
  endpoints: (builder) => ({
    runPipeline: builder.query<
      RunPipelineDTO,
      {
        requestId: string
        elementListForRun: { nodeList: Node[]; edgeList: Edge[] }
        nwbParam: ParamMap
        snakemakeParam: ParamMap
      }
    >({
      // リクエストするたびにキャッシュをクリアするためにidを振っておく
      query: ({ requestId }) => `run/ready/${requestId}`,
      async onCacheEntryAdded(
        { elementListForRun, nwbParam },
        { updateCachedData, cacheDataLoaded, cacheEntryRemoved },
      ) {
        const ws = new WebSocket(`${WS_BASE_URL}/run`)
        try {
          ws.addEventListener('open', () =>
            ws.send(JSON.stringify({ elementListForRun, nwbParam })),
          )
          await cacheDataLoaded
          const listener = (event: MessageEvent) => {
            const data = JSON.parse(event.data)
            updateCachedData((draft) => {
              draft.status = data.status
              draft.message = data.message
              draft.name = data.name
              draft.outputPaths = {
                ...draft.outputPaths,
                ...data.outputPaths,
              }
            })
          }
          ws.addEventListener('message', listener)
        } catch (e) {
          console.log(e)
        }
        await cacheEntryRemoved
        ws.close()
      },
    }),
  }),
})

export const { useLazyRunPipelineQuery, useRunPipelineQuery } = webSocketApi

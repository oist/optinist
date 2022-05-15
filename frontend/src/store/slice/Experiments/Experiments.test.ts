import {
  deleteExperimentByList,
  deleteExperimentByUid,
  getExperiments,
} from './ExperimentsActions'
import reducer, { initialState } from './ExperimentsSlice'

describe('Experiments', () => {
  const uid1 = '96844a59'
  const uid2 = 'eee981f0'
  const getExperimentsPendingAction = {
    type: getExperiments.pending.type,
    meta: {
      requestId: 'FmYmw6sCHA2Ll5JJfPuJN',
      requestStatus: 'pending',
    },
  }
  const getExperimentsFulfilledAction = {
    type: getExperiments.fulfilled.type,
    meta: {
      requestId: 'FmYmw6sCHA2Ll5JJfPuJN',
      requestStatus: 'pending',
    },
    payload: {
      [uid1]: {
        timestamp: '2022-05-07 05:26:54',
        name: 'record test',
        unique_id: uid1,
        function: {
          dummy_image2image8time_4mrz8h7hyk: {
            unique_id: 'dummy_image2image8time_4mrz8h7hyk',
            name: 'dummy_image2image8time',
            success: 'success',
          },
          dummy_image2image_c8tqfxw0mq: {
            unique_id: 'dummy_image2image_c8tqfxw0mq',
            name: 'dummy_image2image',
            success: 'success',
          },
          input_0: {
            unique_id: 'input_0',
            name: 'hoge.tif',
            success: 'success',
          },
        },
        nodeDict: [],
        edgeDict: [],
      },
      [uid2]: {
        timestamp: '2022-05-07 05:54:53',
        name: 'New flow',
        unique_id: uid2,
        function: {
          dummy_image2image8time_4mrz8h7hyk: {
            unique_id: 'dummy_image2image8time_4mrz8h7hyk',
            name: 'dummy_image2image8time',
            success: 'success',
          },
          dummy_image2image_c8tqfxw0mq: {
            unique_id: 'dummy_image2image_c8tqfxw0mq',
            name: 'dummy_image2image',
            success: 'success',
          },
          input_0: {
            unique_id: 'input_0',
            name: 'hoge.tif',
            success: 'success',
          },
        },
        nodeDict: [],
        edgeDict: [],
      },
    },
  }

  test(getExperiments.fulfilled.type, () => {
    const targetState = reducer(
      reducer(initialState, getExperimentsPendingAction),
      getExperimentsFulfilledAction,
    )
    const expectState = {
      status: 'fulfilled',
      experimentList: {
        [uid1]: {
          uid: uid1,
          timestamp: '2022-05-07 05:26:54',
          name: 'record test',
          functions: {
            dummy_image2image8time_4mrz8h7hyk: {
              name: 'dummy_image2image8time',
              nodeId: 'dummy_image2image8time_4mrz8h7hyk',
              status: 'success',
            },
            dummy_image2image_c8tqfxw0mq: {
              name: 'dummy_image2image',
              nodeId: 'dummy_image2image_c8tqfxw0mq',
              status: 'success',
            },
            input_0: { name: 'hoge.tif', nodeId: 'input_0', status: 'success' },
          },
        },
        [uid2]: {
          uid: uid2,
          timestamp: '2022-05-07 05:54:53',
          name: 'New flow',
          functions: {
            dummy_image2image8time_4mrz8h7hyk: {
              name: 'dummy_image2image8time',
              nodeId: 'dummy_image2image8time_4mrz8h7hyk',
              status: 'success',
            },
            dummy_image2image_c8tqfxw0mq: {
              name: 'dummy_image2image',
              nodeId: 'dummy_image2image_c8tqfxw0mq',
              status: 'success',
            },
            input_0: { name: 'hoge.tif', nodeId: 'input_0', status: 'success' },
          },
        },
      },
    }

    expect(targetState).toEqual(expectState)
  })

  // 単一削除
  test(deleteExperimentByUid.fulfilled.type, () => {
    const prevState = reducer(
      reducer(initialState, getExperimentsPendingAction),
      getExperimentsFulfilledAction,
    )
    const deleteExperimentByUidFulfilledAction = {
      type: deleteExperimentByUid.fulfilled.type,
      payload: true,
      meta: {
        arg: uid2,
        requestId: 'FZcujcGG38JPwbQDR6c-J',
        requestStatus: 'fulfilled',
      },
    }
    const targetState = reducer(prevState, deleteExperimentByUidFulfilledAction)
    expect(
      targetState.status === 'fulfilled' &&
        !Object.keys(targetState.experimentList).includes(uid2),
    ).toBe(true)
  })

  // 複数削除
  test(deleteExperimentByList.fulfilled.type, () => {
    const prevState = reducer(
      reducer(initialState, getExperimentsPendingAction),
      getExperimentsFulfilledAction,
    )
    const deleteExperimentByListFulfilledAction = {
      type: deleteExperimentByList.fulfilled.type,
      payload: true,
      meta: {
        arg: [uid1, uid2],
        requestId: 'faNrL5ZV3SRODugFlJM20',
        requestStatus: 'fulfilled',
      },
    }
    const targetState = reducer(
      prevState,
      deleteExperimentByListFulfilledAction,
    )
    expect(
      targetState.status === 'fulfilled' && targetState.experimentList,
    ).toEqual({})
  })
})

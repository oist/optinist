import { getAlgoList } from './AlgorithmListActions'
import reducer, { initialState } from './AlgorithmListSlice'
import { AlgorithmListType } from './AlgorithmListType'

describe('getAlgoList', () => {
  test(getAlgoList.fulfilled.type, () => {
    expect(
      reducer(initialState, {
        type: getAlgoList.fulfilled.type,
        payload: mockPayload,
        meta: {
          requestId: 'jch5KyN83_x3ZtGQssSnJ',
          requestStatus: 'fulfilled',
        },
      }),
    ).toEqual(expectState)
  })
  const mockPayload = {
    caiman: {
      children: {
        caiman_mc: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'mc_images', type: 'ImageData' }],
          parameter: 'caiman_mc.yaml',
          path: 'caiman/caiman_mc',
        },
        caiman_cnmf: {
          args: [{ name: 'images', type: 'ImageData', isNone: false }],
          returns: [
            { name: 'fluorescence', type: 'FluoData' },
            { name: 'iscell', type: 'IscellData' },
          ],
          parameter: 'caiman_cnmf.yaml',
          path: 'caiman/caiman_cnmf',
        },
      },
    },
    suite2p: {
      children: {
        suite2p_file_convert: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'ops', type: 'Suite2pData' }],
          parameter: null,
          path: 'suite2p/suite2p_file_convert',
        },
        suite2p_registration: {
          args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
          returns: [{ name: 'ops', type: 'Suite2pData' }],
          parameter: null,
          path: 'suite2p/suite2p_registration',
        },
        suite2p_roi: {
          args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
          returns: [
            { name: 'ops', type: 'Suite2pData' },
            { name: 'fluorescence', type: 'FluoData' },
            { name: 'iscell', type: 'IscellData' },
          ],
          parameter: null,
          path: 'suite2p/suite2p_roi',
        },
        suite2p_spike_deconv: {
          args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
          returns: [
            { name: 'ops', type: 'Suite2pData' },
            { name: 'spks', type: 'FluoData' },
          ],
          parameter: null,
          path: 'suite2p/suite2p_spike_deconv',
        },
      },
    },
    dummy: {
      children: {
        dummy_image2image: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'image2image', type: 'ImageData' }],
          parameter: null,
          path: 'dummy/dummy_image2image',
        },
        dummy_image2time: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'image2time', type: 'TimeSeriesData' }],
          parameter: null,
          path: 'dummy/dummy_image2time',
        },
        dummy_image2heat: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'image2heat', type: 'HeatMapData' }],
          parameter: null,
          path: 'dummy/dummy_image2heat',
        },
        dummy_time2time: {
          args: [{ name: 'timeseries', type: 'TimeSeriesData', isNone: false }],
          returns: [{ name: 'time2time', type: 'TimeSeriesData' }],
          parameter: null,
          path: 'dummy/dummy_time2time',
        },
        dummy_image2image8time: {
          args: [{ name: 'image1', type: 'ImageData', isNone: false }],
          returns: [
            { name: 'image', type: 'ImageData' },
            { name: 'timeseries', type: 'TimeSeriesData' },
          ],
          parameter: null,
          path: 'dummy/dummy_image2image8time',
        },
        dummy_keyerror: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [],
          parameter: null,
          path: 'dummy/dummy_keyerror',
        },
        dummy_typeerror: {
          args: [{ name: 'image', type: 'str', isNone: false }],
          returns: [],
          parameter: null,
          path: 'dummy/dummy_typeerror',
        },
        dummy_image2time8iscell: {
          args: [{ name: 'image1', type: 'ImageData', isNone: false }],
          returns: [
            { name: 'timeseries', type: 'TimeSeriesData' },
            { name: 'iscell', type: 'IscellData' },
          ],
          parameter: null,
          path: 'dummy/dummy_image2time8iscell',
        },
        dummy_image2roi: {
          args: [{ name: 'image1', type: 'ImageData', isNone: false }],
          returns: [{ name: 'roi', type: 'RoiData' }],
          parameter: null,
          path: 'dummy/dummy_image2roi',
        },
        dummy_image2image8roi: {
          args: [{ name: 'image1', type: 'ImageData', isNone: false }],
          returns: [
            { name: 'image', type: 'ImageData' },
            { name: 'roi', type: 'RoiData' },
          ],
          parameter: null,
          path: 'dummy/dummy_image2image8roi',
        },
        dummy_image2image8roi8time8heat: {
          args: [{ name: 'image1', type: 'ImageData', isNone: false }],
          returns: [
            { name: 'image', type: 'ImageData' },
            { name: 'roi', type: 'RoiData' },
            { name: 'timeseries', type: 'TimeSeriesData' },
            { name: 'heat', type: 'HeatMapData' },
          ],
          parameter: null,
          path: 'dummy/dummy_image2image8roi8time8heat',
        },
        dummy_image2scatter: {
          args: [{ name: 'image', type: 'ImageData', isNone: false }],
          returns: [{ name: 'scatter', type: 'ScatterData' }],
          parameter: null,
          path: 'dummy/dummy_image2scatter',
        },
      },
    },
    studio: {
      children: {
        basic_neural_analysis: {
          children: {
            eta: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                {
                  name: 'behaviors_data',
                  type: 'BehaviorData',
                  isNone: false,
                },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [{ name: 'mean', type: 'TimeSeriesData' }],
              parameter: null,
              path: 'studio/basic_neural_analysis/eta',
            },
            cell_grouping: {
              args: [
                {
                  name: 'neural_data',
                  type: 'TimeSeriesData',
                  isNone: false,
                },
              ],
              returns: [],
              parameter: null,
              path: 'studio/basic_neural_analysis/cell_grouping',
            },
          },
        },
        dimension_reduction: {
          children: {
            cca: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                {
                  name: 'behaviors_data',
                  type: 'BehaviorData',
                  isNone: false,
                },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/dimension_reduction/cca',
            },
            pca: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/dimension_reduction/pca',
            },
            tsne: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/dimension_reduction/tsne',
            },
          },
        },
        neural_population_analysis: {
          children: {
            correlation: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_population_analysis/correlation',
            },
            cross_correlation: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_population_analysis/cross_correlation',
            },
            granger: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_population_analysis/granger',
            },
          },
        },
        neural_decoding: {
          children: {
            glm: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                {
                  name: 'behaviors_data',
                  type: 'BehaviorData',
                  isNone: false,
                },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_decoding/glm',
            },
            lda: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                {
                  name: 'behaviors_data',
                  type: 'BehaviorData',
                  isNone: false,
                },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_decoding/lda',
            },
            svm: {
              args: [
                { name: 'neural_data', type: 'FluoData', isNone: false },
                {
                  name: 'behaviors_data',
                  type: 'BehaviorData',
                  isNone: false,
                },
                { name: 'iscell', type: 'IscellData', isNone: true },
              ],
              returns: [],
              parameter: null,
              path: 'studio/neural_decoding/svm',
            },
          },
        },
      },
    },
  }
  const expectState: AlgorithmListType = {
    isLatest: true,
    tree: {
      caiman: {
        type: 'parent',
        children: {
          caiman_mc: {
            type: 'child',
            functionPath: 'caiman/caiman_mc',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'mc_images', type: 'ImageData' }],
          },
          caiman_cnmf: {
            type: 'child',
            functionPath: 'caiman/caiman_cnmf',
            args: [{ name: 'images', type: 'ImageData', isNone: false }],
            returns: [
              { name: 'fluorescence', type: 'FluoData' },
              { name: 'iscell', type: 'IscellData' },
            ],
          },
        },
      },
      suite2p: {
        type: 'parent',
        children: {
          suite2p_file_convert: {
            type: 'child',
            functionPath: 'suite2p/suite2p_file_convert',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'ops', type: 'Suite2pData' }],
          },
          suite2p_registration: {
            type: 'child',
            functionPath: 'suite2p/suite2p_registration',
            args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
            returns: [{ name: 'ops', type: 'Suite2pData' }],
          },
          suite2p_roi: {
            type: 'child',
            functionPath: 'suite2p/suite2p_roi',
            args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
            returns: [
              { name: 'ops', type: 'Suite2pData' },
              { name: 'fluorescence', type: 'FluoData' },
              { name: 'iscell', type: 'IscellData' },
            ],
          },
          suite2p_spike_deconv: {
            type: 'child',
            functionPath: 'suite2p/suite2p_spike_deconv',
            args: [{ name: 'ops', type: 'Suite2pData', isNone: false }],
            returns: [
              { name: 'ops', type: 'Suite2pData' },
              { name: 'spks', type: 'FluoData' },
            ],
          },
        },
      },
      dummy: {
        type: 'parent',
        children: {
          dummy_image2image: {
            type: 'child',
            functionPath: 'dummy/dummy_image2image',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'image2image', type: 'ImageData' }],
          },
          dummy_image2time: {
            type: 'child',
            functionPath: 'dummy/dummy_image2time',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'image2time', type: 'TimeSeriesData' }],
          },
          dummy_image2heat: {
            type: 'child',
            functionPath: 'dummy/dummy_image2heat',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'image2heat', type: 'HeatMapData' }],
          },
          dummy_time2time: {
            type: 'child',
            functionPath: 'dummy/dummy_time2time',
            args: [
              { name: 'timeseries', type: 'TimeSeriesData', isNone: false },
            ],
            returns: [{ name: 'time2time', type: 'TimeSeriesData' }],
          },
          dummy_image2image8time: {
            type: 'child',
            functionPath: 'dummy/dummy_image2image8time',
            args: [{ name: 'image1', type: 'ImageData', isNone: false }],
            returns: [
              { name: 'image', type: 'ImageData' },
              { name: 'timeseries', type: 'TimeSeriesData' },
            ],
          },
          dummy_keyerror: {
            type: 'child',
            functionPath: 'dummy/dummy_keyerror',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [],
          },
          dummy_typeerror: {
            type: 'child',
            functionPath: 'dummy/dummy_typeerror',
            args: [{ name: 'image', type: 'str', isNone: false }],
            returns: [],
          },
          dummy_image2time8iscell: {
            type: 'child',
            functionPath: 'dummy/dummy_image2time8iscell',
            args: [{ name: 'image1', type: 'ImageData', isNone: false }],
            returns: [
              { name: 'timeseries', type: 'TimeSeriesData' },
              { name: 'iscell', type: 'IscellData' },
            ],
          },
          dummy_image2roi: {
            type: 'child',
            functionPath: 'dummy/dummy_image2roi',
            args: [{ name: 'image1', type: 'ImageData', isNone: false }],
            returns: [{ name: 'roi', type: 'RoiData' }],
          },
          dummy_image2image8roi: {
            type: 'child',
            functionPath: 'dummy/dummy_image2image8roi',
            args: [{ name: 'image1', type: 'ImageData', isNone: false }],
            returns: [
              { name: 'image', type: 'ImageData' },
              { name: 'roi', type: 'RoiData' },
            ],
          },
          dummy_image2image8roi8time8heat: {
            type: 'child',
            functionPath: 'dummy/dummy_image2image8roi8time8heat',
            args: [{ name: 'image1', type: 'ImageData', isNone: false }],
            returns: [
              { name: 'image', type: 'ImageData' },
              { name: 'roi', type: 'RoiData' },
              { name: 'timeseries', type: 'TimeSeriesData' },
              { name: 'heat', type: 'HeatMapData' },
            ],
          },
          dummy_image2scatter: {
            type: 'child',
            functionPath: 'dummy/dummy_image2scatter',
            args: [{ name: 'image', type: 'ImageData', isNone: false }],
            returns: [{ name: 'scatter', type: 'ScatterData' }],
          },
        },
      },
      studio: {
        type: 'parent',
        children: {
          basic_neural_analysis: {
            type: 'parent',
            children: {
              eta: {
                type: 'child',
                functionPath: 'studio/basic_neural_analysis/eta',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  {
                    name: 'behaviors_data',
                    type: 'BehaviorData',
                    isNone: false,
                  },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [{ name: 'mean', type: 'TimeSeriesData' }],
              },
              cell_grouping: {
                type: 'child',
                functionPath: 'studio/basic_neural_analysis/cell_grouping',
                args: [
                  {
                    name: 'neural_data',
                    type: 'TimeSeriesData',
                    isNone: false,
                  },
                ],
                returns: [],
              },
            },
          },
          dimension_reduction: {
            type: 'parent',
            children: {
              cca: {
                type: 'child',
                functionPath: 'studio/dimension_reduction/cca',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  {
                    name: 'behaviors_data',
                    type: 'BehaviorData',
                    isNone: false,
                  },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              pca: {
                type: 'child',
                functionPath: 'studio/dimension_reduction/pca',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              tsne: {
                type: 'child',
                functionPath: 'studio/dimension_reduction/tsne',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
            },
          },
          neural_population_analysis: {
            type: 'parent',
            children: {
              correlation: {
                type: 'child',
                functionPath: 'studio/neural_population_analysis/correlation',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              cross_correlation: {
                type: 'child',
                functionPath:
                  'studio/neural_population_analysis/cross_correlation',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              granger: {
                type: 'child',
                functionPath: 'studio/neural_population_analysis/granger',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
            },
          },
          neural_decoding: {
            type: 'parent',
            children: {
              glm: {
                type: 'child',
                functionPath: 'studio/neural_decoding/glm',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  {
                    name: 'behaviors_data',
                    type: 'BehaviorData',
                    isNone: false,
                  },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              lda: {
                type: 'child',
                functionPath: 'studio/neural_decoding/lda',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  {
                    name: 'behaviors_data',
                    type: 'BehaviorData',
                    isNone: false,
                  },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
              svm: {
                type: 'child',
                functionPath: 'studio/neural_decoding/svm',
                args: [
                  { name: 'neural_data', type: 'FluoData', isNone: false },
                  {
                    name: 'behaviors_data',
                    type: 'BehaviorData',
                    isNone: false,
                  },
                  { name: 'iscell', type: 'IscellData', isNone: true },
                ],
                returns: [],
              },
            },
          },
        },
      },
    },
  }
})

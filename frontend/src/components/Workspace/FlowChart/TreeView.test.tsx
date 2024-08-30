/* eslint-disable no-undef */
import React from "react"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"
import thunk from "redux-thunk"

import { describe, it, beforeEach } from "@jest/globals"
import { render, screen, fireEvent } from "@testing-library/react"
import { userEvent } from "@testing-library/user-event"

import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"

jest.mock("store/slice/AlgorithmList/AlgorithmListActions", () => ({
  getAlgoList: jest.fn(),
}))

jest.mock("react-dnd", () => ({
  ...jest.requireActual("react-dnd"),
  useDrag: () => [{ isDragging: false }, null],
}))

const mockStore = configureStore([thunk])

describe("AlgorithmTreeView", () => {
  let store: ReturnType<typeof mockStore>

  beforeEach(() => {
    store = mockStore({
      isLatest: true,
      algorithmList: {
        tree: {
          caiman: {
            type: "parent",
            children: {
              caiman_mc: {
                type: "child",
                functionPath: "caiman/caiman_mc",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "mc_images", type: "ImageData" }],
              },
              caiman_cnmf: {
                type: "child",
                functionPath: "caiman/caiman_cnmf",
                args: [{ name: "images", type: "ImageData", isNone: false }],
                returns: [
                  { name: "fluorescence", type: "FluoData" },
                  { name: "iscell", type: "IscellData" },
                ],
              },
            },
          },
          suite2p: {
            type: "parent",
            children: {
              suite2p_file_convert: {
                type: "child",
                functionPath: "suite2p/suite2p_file_convert",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "ops", type: "Suite2pData" }],
              },
              suite2p_registration: {
                type: "child",
                functionPath: "suite2p/suite2p_registration",
                args: [{ name: "ops", type: "Suite2pData", isNone: false }],
                returns: [{ name: "ops", type: "Suite2pData" }],
              },
              suite2p_roi: {
                type: "child",
                functionPath: "suite2p/suite2p_roi",
                args: [{ name: "ops", type: "Suite2pData", isNone: false }],
                returns: [
                  { name: "ops", type: "Suite2pData" },
                  { name: "fluorescence", type: "FluoData" },
                  { name: "iscell", type: "IscellData" },
                ],
              },
              suite2p_spike_deconv: {
                type: "child",
                functionPath: "suite2p/suite2p_spike_deconv",
                args: [{ name: "ops", type: "Suite2pData", isNone: false }],
                returns: [
                  { name: "ops", type: "Suite2pData" },
                  { name: "spks", type: "FluoData" },
                ],
              },
            },
          },
          lccd: {
            type: "parent",
            children: {
              lccd_image2image: {
                type: "child",
                functionPath: "lccd/lccd_image2image",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "image2image", type: "ImageData" }],
              },
              lccd_image2time: {
                type: "child",
                functionPath: "lccd/lccd_image2time",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "image2time", type: "TimeSeriesData" }],
              },
              lccd_image2heat: {
                type: "child",
                functionPath: "lccd/lccd_image2heat",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "image2heat", type: "HeatMapData" }],
              },
              lccd_time2time: {
                type: "child",
                functionPath: "lccd/lccd_time2time",
                args: [
                  { name: "timeseries", type: "TimeSeriesData", isNone: false },
                ],
                returns: [{ name: "time2time", type: "TimeSeriesData" }],
              },
              lccd_image2image8time: {
                type: "child",
                functionPath: "lccd/lccd_image2image8time",
                args: [{ name: "image1", type: "ImageData", isNone: false }],
                returns: [
                  { name: "image", type: "ImageData" },
                  { name: "timeseries", type: "TimeSeriesData" },
                ],
              },
              lccd_keyerror: {
                type: "child",
                functionPath: "lccd/lccd_keyerror",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [],
              },
              lccd_typeerror: {
                type: "child",
                functionPath: "lccd/lccd_typeerror",
                args: [{ name: "image", type: "str", isNone: false }],
                returns: [],
              },
              lccd_image2time8iscell: {
                type: "child",
                functionPath: "lccd/lccd_image2time8iscell",
                args: [{ name: "image1", type: "ImageData", isNone: false }],
                returns: [
                  { name: "timeseries", type: "TimeSeriesData" },
                  { name: "iscell", type: "IscellData" },
                ],
              },
              lccd_image2roi: {
                type: "child",
                functionPath: "lccd/lccd_image2roi",
                args: [{ name: "image1", type: "ImageData", isNone: false }],
                returns: [{ name: "roi", type: "RoiData" }],
              },
              lccd_image2image8roi: {
                type: "child",
                functionPath: "lccd/lccd_image2image8roi",
                args: [{ name: "image1", type: "ImageData", isNone: false }],
                returns: [
                  { name: "image", type: "ImageData" },
                  { name: "roi", type: "RoiData" },
                ],
              },
              lccd_image2image8roi8time8heat: {
                type: "child",
                functionPath: "lccd/lccd_image2image8roi8time8heat",
                args: [{ name: "image1", type: "ImageData", isNone: false }],
                returns: [
                  { name: "image", type: "ImageData" },
                  { name: "roi", type: "RoiData" },
                  { name: "timeseries", type: "TimeSeriesData" },
                  { name: "heat", type: "HeatMapData" },
                ],
              },
              lccd_image2scatter: {
                type: "child",
                functionPath: "lccd/lccd_image2scatter",
                args: [{ name: "image", type: "ImageData", isNone: false }],
                returns: [{ name: "scatter", type: "ScatterData" }],
              },
            },
          },
          optinist: {
            type: "parent",
            children: {
              basic_neural_analysis: {
                type: "parent",
                children: {
                  eta: {
                    type: "child",
                    functionPath: "optinist/basic_neural_analysis/eta",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      {
                        name: "behaviors_data",
                        type: "BehaviorData",
                        isNone: false,
                      },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [{ name: "mean", type: "TimeSeriesData" }],
                  },
                  cell_grouping: {
                    type: "child",
                    functionPath:
                      "optinist/basic_neural_analysis/cell_grouping",
                    args: [
                      {
                        name: "neural_data",
                        type: "TimeSeriesData",
                        isNone: false,
                      },
                    ],
                    returns: [],
                  },
                },
              },
              dimension_reduction: {
                type: "parent",
                children: {
                  cca: {
                    type: "child",
                    functionPath: "optinist/dimension_reduction/cca",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      {
                        name: "behaviors_data",
                        type: "BehaviorData",
                        isNone: false,
                      },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  pca: {
                    type: "child",
                    functionPath: "optinist/dimension_reduction/pca",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  tsne: {
                    type: "child",
                    functionPath: "optinist/dimension_reduction/tsne",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                },
              },
              neural_population_analysis: {
                type: "parent",
                children: {
                  correlation: {
                    type: "child",
                    functionPath:
                      "optinist/neural_population_analysis/correlation",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  cross_correlation: {
                    type: "child",
                    functionPath:
                      "optinist/neural_population_analysis/cross_correlation",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  granger: {
                    type: "child",
                    functionPath: "optinist/neural_population_analysis/granger",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                },
              },
              neural_decoding: {
                type: "parent",
                children: {
                  glm: {
                    type: "child",
                    functionPath: "optinist/neural_decoding/glm",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      {
                        name: "behaviors_data",
                        type: "BehaviorData",
                        isNone: false,
                      },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  lda: {
                    type: "child",
                    functionPath: "optinist/neural_decoding/lda",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      {
                        name: "behaviors_data",
                        type: "BehaviorData",
                        isNone: false,
                      },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                  svm: {
                    type: "child",
                    functionPath: "optinist/neural_decoding/svm",
                    args: [
                      { name: "neural_data", type: "FluoData", isNone: false },
                      {
                        name: "behaviors_data",
                        type: "BehaviorData",
                        isNone: false,
                      },
                      { name: "iscell", type: "IscellData", isNone: true },
                    ],
                    returns: [],
                  },
                },
              },
            },
          },
        },
      },
      pipeline: {
        currentPipeline: {
          uid: "123",
          name: "pipeline1",
        },
      },
    })
    store.dispatch = jest.fn()
  })

  it("renders the AlgorithmTreeView component", async () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )
    expect(screen.getByText("Data")).toBeInTheDocument()
    expect(screen.getByText("Algorithm")).toBeInTheDocument()
  })

  it("renders the AlgorithmTree Data TreeItems", async () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    // Click on the "Data" TreeItem to expand it
    const dataTreeLabel = screen.getByText("Data")
    await userEvent.click(dataTreeLabel)

    // Check if all the TreeItems are rendered
    expect(screen.getByText("image")).toBeInTheDocument()
    expect(screen.getByText("csv")).toBeInTheDocument()
    expect(screen.getByText("hdf5")).toBeInTheDocument()
    expect(screen.getByText("fluo")).toBeInTheDocument()
    expect(screen.getByText("behavior")).toBeInTheDocument()
    expect(screen.getByText("matlab")).toBeInTheDocument()
    expect(screen.getByText("microscope")).toBeInTheDocument()
  })

  it("renders the AlgorithmTree Algorithm TreeItems", async () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    // Click on the "Data" TreeItem to expand it
    const algorithmTreeLabel = screen.getByText("Algorithm")
    await userEvent.click(algorithmTreeLabel)

    // Check if all the TreeItems are rendered
    expect(screen.getByText("caiman")).toBeInTheDocument()
    expect(screen.getByText("suite2p")).toBeInTheDocument()
    expect(screen.getByText("lccd")).toBeInTheDocument()
    expect(screen.getByText("optinist")).toBeInTheDocument()
  })

  it("dispatches getAlgoList action when component mounts", () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    expect(store.dispatch).toHaveBeenCalledWith(getAlgoList())
  })

  it("dispatches the correct action when the add button is clicked", async () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    const dataTreeLabel = screen.getByText("Data")
    await userEvent.click(dataTreeLabel)

    // Click on the "image" TreeItem to select it
    const addButton = screen.getAllByLabelText("add")[0]
    await userEvent.click(addButton)

    // Verify that the action was dispatched with the expected payload
    expect(store.dispatch).toHaveBeenCalledWith(
      expect.objectContaining({
        type: "flowElement/addInputNode",
        payload: {
          node: expect.any(Object),
          fileType: "image",
        },
      }),
    )
  })
})

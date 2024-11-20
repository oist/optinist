/* eslint-disable no-undef */
import React from "react"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"
import thunk from "redux-thunk"

import { describe, it, beforeEach } from "@jest/globals"
import { prettyDOM, render, screen } from "@testing-library/react"
import { userEvent } from "@testing-library/user-event"

import { mockStoreData } from "components/Workspace/__tests__/FlowChart/mockStoreData.json"
import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"
import { addAlgorithmNode } from "store/slice/FlowElement/FlowElementActions"

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
    store = mockStore(mockStoreData)
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

  it("dispatches the correct action when the Image node add button is clicked", async () => {
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

  it("dispatches the correct action when the algorithm add button is clicked", async () => {
    // Spy on the `addAlgorithmNode` action for this test
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    const addAlgorithmNodeMock = jest
      .spyOn(
        // eslint-disable-next-line @typescript-eslint/no-var-requires
        require("store/slice/FlowElement/FlowElementActions"),
        "addAlgorithmNode",
      )
      .mockImplementation(jest.fn())

    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    // Click the "Algorithm" label to expand the node
    const algorithmTreeLabel = screen.getByText("Algorithm")
    await userEvent.click(algorithmTreeLabel)

    // Ensure the "caiman" node exists and click it
    const caimanTreeLabel = screen.getByText("caiman")
    await userEvent.click(caimanTreeLabel)

    // Click the add button for the "caiman" node (using a more specific selector if needed)
    const addButton = screen.getAllByTestId("AddIcon")[0]
    await userEvent.click(addButton)

    // Verify that the addAlgorithmNode action was dispatched
    expect(addAlgorithmNode).toHaveBeenCalledWith({
      node: {
        data: { label: "caiman_mc", type: "algorithm" },
        id: expect.any(String), // Use `expect.any(String)` if the id is dynamically generated,
        position: undefined,
        type: "AlgorithmNode",
      },
      name: "caiman_mc",
      functionPath: "caiman/caiman_mc",
      runAlready: true,
    })
  })
})

/* eslint-disable no-undef */
import React from "react"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"
import thunk from "redux-thunk"

import { describe, it, beforeEach } from "@jest/globals"
import { render, screen } from "@testing-library/react"
import { userEvent } from "@testing-library/user-event"

import { mockStoreData } from "components/Workspace/FlowChart/testdata/mockStoreData"
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

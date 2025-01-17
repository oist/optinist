import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { Store, AnyAction } from "@reduxjs/toolkit"
import { render, screen, fireEvent } from "@testing-library/react"

import { TreeItemLabel } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import { deleteFile } from "store/slice/FilesTree/FilesTreeAction"
import { AppDispatch } from "store/store"

jest.mock("store/slice/FilesTree/FilesTreeAction", () => ({
  deleteFile: jest.fn(),
}))

const mockStore = configureStore<
  Partial<{ workspace: { currentWorkspace: { workspaceId?: number } } }>,
  AppDispatch
>([])

describe("TreeItemLabel Component", () => {
  let store: Store<unknown, AnyAction>

  beforeEach(() => {
    store = mockStore({
      workspace: {
        currentWorkspace: {
          workspaceId: 123,
        },
      },
    })
    store.dispatch = jest.fn()
  })

  it("should dispatch deleteFile action when delete is confirmed", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: false, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    // Simulate opening the delete confirmation dialog
    const deleteIcon = screen.getByTestId("DeleteIcon")
    fireEvent.click(deleteIcon)

    // Confirm the delete action
    const confirmButton = screen.getByRole("button", { name: /DELETE/i })
    fireEvent.click(confirmButton)

    expect(store.dispatch).toHaveBeenCalledWith(
      deleteFile({
        workspaceId: 123,
        fileName: "testFile",
        fileType: "image",
      }),
    )
  })

  it("should not dispatch deleteFile action if confirmation clicked no", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: false, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    // Simulate opening the delete confirmation dialog
    const deleteIcon = screen.getByTestId("DeleteIcon")
    fireEvent.click(deleteIcon)

    // Confirm the delete action
    const confirmButton = screen.getByRole("button", { name: /CANCEL/i })
    fireEvent.click(confirmButton)

    expect(store.dispatch).not.toHaveBeenCalledWith(deleteFile)
  })

  it("should disable delete button if the file checkbox is checked", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: true, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    const deleteButton = screen.getByTestId("DeleteIconBtn")
    expect(deleteButton).toBeDisabled()
  })

  it("should enable delete button if the file checkbox is not checked", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: false, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    const deleteButton = screen.getByTestId("DeleteIcon")
    expect(deleteButton).toBeEnabled()
  })
})

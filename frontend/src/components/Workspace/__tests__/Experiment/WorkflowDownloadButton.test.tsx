/* eslint-disable no-undef */
import "@testing-library/jest-dom"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { render, screen, fireEvent, waitFor } from "@testing-library/react"

import { downloadWorkflowConfigApi } from "api/workflow/Workflow"
import mockStoreData from "components/Workspace/__tests__/Experiment/mockStore.json"
import { WorkflowDownloadButton } from "components/Workspace/Experiment/Button/DownloadButton"
import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable"

jest.mock("api/workflow/Workflow", () => ({
  downloadWorkflowConfigApi: jest.fn(),
}))

jest.mock("notistack", () => ({
  useSnackbar: () => ({
    enqueueSnackbar: jest.fn(),
  }),
}))

const mockStore = configureStore([])
const store = mockStore(mockStoreData)

describe("WorkflowDownloadButton", () => {
  it("triggers file download when the workflow download button is clicked", async () => {
    render(
      <Provider store={store}>
        <ExperimentUidContext.Provider value="exp1">
          <WorkflowDownloadButton />
        </ExperimentUidContext.Provider>
      </Provider>,
    )

    // Find the workflow download button using data-testid
    const downloadButton = screen.getByTestId("workflow-download-button")

    // Mock the click event
    fireEvent.click(downloadButton)

    // Wait for the download URL to be created and set
    await waitFor(() => {
      expect(downloadWorkflowConfigApi).toHaveBeenCalledWith(1, "exp1")
    })

    // Get the hidden anchor element used for downloading using data-testid
    const link = await waitFor(() =>
      screen.getByTestId("workflow-download-link"),
    )

    // Check that the download attribute is set correctly
    expect(link).toHaveAttribute("download", "workflow_exp1.yaml")
  })
})

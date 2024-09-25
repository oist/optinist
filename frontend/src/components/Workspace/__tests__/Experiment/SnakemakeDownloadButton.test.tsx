import "@testing-library/jest-dom"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { render, screen, fireEvent, waitFor } from "@testing-library/react"

import { downloadExperimentConfigApi } from "api/experiments/Experiments"
import mockStoreData from "components/Workspace/__tests__/Experiment/mockStore.json"
import { SnakemakeDownloadButton } from "components/Workspace/Experiment/Button/DownloadButton"
import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable"

// Mock the API module
jest.mock("api/experiments/Experiments", () => ({
  downloadExperimentConfigApi: jest.fn(),
}))

jest.mock("notistack", () => ({
  useSnackbar: () => ({
    enqueueSnackbar: jest.fn(),
  }),
}))

const mockStore = configureStore([])
const store = mockStore(mockStoreData)

describe("SnakemakeDownloadButton", () => {
  beforeEach(() => {
    // Mock implementation of the API function
    ;(downloadExperimentConfigApi as jest.Mock).mockResolvedValue(
      new Blob(["mocked file content"], { type: "application/yaml" }),
    )
  })

  it("checks if the Snakemake download button is clickable and triggers file download", async () => {
    render(
      <Provider store={store}>
        <ExperimentUidContext.Provider value="exp1">
          <SnakemakeDownloadButton />
        </ExperimentUidContext.Provider>
      </Provider>,
    )

    // Find the Snakemake download button
    const downloadButton = screen.getByRole("button")

    // Check if the button is enabled (i.e., clickable)
    expect(downloadButton).toBeEnabled()

    // Mock the click event
    fireEvent.click(downloadButton)

    // Wait for the downloadExperimentConfigApi to be called
    await waitFor(() => {
      expect(downloadExperimentConfigApi).toHaveBeenCalledWith(1, "exp1")
    })

    // Get the hidden anchor element using data-testid
    const link = await waitFor(() =>
      screen.getByTestId("snakemake-download-link"),
    )

    // Check that the download attribute is set correctly
    expect(link).toHaveAttribute("download", "snakemake_exp1.yaml")
  })
})

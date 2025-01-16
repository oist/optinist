import "@testing-library/jest-dom"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { render, screen, fireEvent, waitFor } from "@testing-library/react"

import { downloadExperimentNwbApi } from "api/experiments/Experiments" // Adjust import path
import mockStoreData from "components/Workspace/__tests__/Experiment/mockStore.json"
import { NWBDownloadButton } from "components/Workspace/Experiment/Button/DownloadButton" // Adjust import path
import { ExperimentUidContext } from "components/Workspace/Experiment/ExperimentTable" // Adjust import path

jest.mock("api/experiments/Experiments", () => ({
  downloadExperimentNwbApi: jest.fn(),
}))

jest.mock("notistack", () => ({
  useSnackbar: () => ({
    enqueueSnackbar: jest.fn(),
  }),
}))

const mockStore = configureStore([])
const store = mockStore(mockStoreData)

describe("NWBDownloadButton", () => {
  beforeEach(() => {
    ;(downloadExperimentNwbApi as jest.Mock).mockResolvedValue(
      "mocked NWB file content",
    )
  })

  it("triggers NWB file download when the button is clicked and hasNWB is true", async () => {
    render(
      <Provider store={store}>
        <ExperimentUidContext.Provider value="exp1">
          <NWBDownloadButton name="testName" nodeId="node1" hasNWB={true} />
        </ExperimentUidContext.Provider>
      </Provider>,
    )

    // Find the NWB download button
    const downloadButton = screen.getByRole("button")

    // Check if the button is enabled (i.e., clickable)
    expect(downloadButton).toBeEnabled()

    // Mock the click event
    fireEvent.click(downloadButton)

    // Wait for the downloadExperimentNwbApi to be called
    await waitFor(() => {
      expect(downloadExperimentNwbApi).toHaveBeenCalledWith(1, "exp1", "node1")
    })

    const link = await waitFor(() => screen.getByTestId("nwb-download-link"))

    // Check that the download attribute is set correctly
    expect(link).toHaveAttribute("download", "nwb_testName.nwb")
  })
})

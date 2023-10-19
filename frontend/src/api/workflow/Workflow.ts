import { ExperimentDTO } from "api/experiments/Experiments"
import { EdgeDict, NodeDict } from "api/run/Run"
import { BASE_URL } from "const/API"
import axios from "utils/axios"


export type WorkflowConfigDTO = {
  nodeDict: NodeDict
  edgeDict: EdgeDict
}
export type WorkflowWithResultDTO = ExperimentDTO & WorkflowConfigDTO

export async function fetchWorkflowApi(
  workspace_id: number,
): Promise<WorkflowWithResultDTO> {
  const response = await axios.get(`${BASE_URL}/workflow/fetch/${workspace_id}`)
  return response.data
}

export async function reproduceWorkflowApi(
  workspaceId: number,
  uid: string,
): Promise<WorkflowWithResultDTO> {
  const response = await axios.get(
    `${BASE_URL}/workflow/reproduce/${workspaceId}/${uid}`,
  )
  return response.data
}

export async function downloadWorkflowConfigApi(
  workspaceId: number,
  uid: string,
) {
  const response = await axios.get(
    `${BASE_URL}/workflow/download/${workspaceId}/${uid}`,
    {
      responseType: "blob",
    },
  )
  return response.data
}

export async function importWorkflowConfigApi(
  formData: FormData,
): Promise<WorkflowConfigDTO> {
  const response = await axios.post(`${BASE_URL}/workflow/import`, formData)
  return response.data
}

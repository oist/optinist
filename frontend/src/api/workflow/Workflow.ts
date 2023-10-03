import axios from 'utils/axios'

import { EdgeDict, NodeDict } from 'api/run/Run'
import { BASE_URL } from 'const/API'

export type WorkflowConfigDTO = {
  nodeDict: NodeDict
  edgeDict: EdgeDict
}

export async function reproduceWorkflowApi(
  workspaceId: number,
  uid: string,
): Promise<WorkflowConfigDTO> {
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
      responseType: 'blob',
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

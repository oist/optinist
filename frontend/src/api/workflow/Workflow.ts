import axios from 'utils/axios'

import { RunPostData } from 'api/run/Run'
import { BASE_URL } from 'const/API'

export async function importWorkflowByUidApi(
  workspaceId: string,
  uid: string,
): Promise<RunPostData> {
  const response = await axios.get(
    `${BASE_URL}/workflow/import/${workspaceId}/${uid}`,
  )
  return response.data
}

import { Connection } from "reactflow"

export function toHandleId(nodeId: string, name: string, type: string) {
  return `${nodeId}--${name}--${type}`
}

export function getHandleType(handleId: string) {
  return handleId.split("--")[2]
}

export function isValidConnection(connection: Connection) {
  if (connection.sourceHandle != null && connection.targetHandle != null) {
    const source = getHandleType(connection.sourceHandle)
    const target = getHandleType(connection.targetHandle)
    // NOTE: SpikingActivityData is the same as FluoData. Just renamed for suite2p_spike_deconv.
    if (source === "SpikingActivityData") {
      return target === "FluoData" || target === "SpikingActivityData"
    } else {
      return source === target
    }
  } else {
    return true
  }
}

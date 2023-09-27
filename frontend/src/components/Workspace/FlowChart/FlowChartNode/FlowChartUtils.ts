import { Connection } from 'react-flow-renderer'

export function toHandleId(nodeId: string, name: string, type: string) {
  return `${nodeId}--${name}--${type}`
}

export function getHandleType(handleId: string) {
  return handleId.split('--')[2]
}

export function isValidConnection(connection: Connection) {
  if (connection.sourceHandle != null && connection.targetHandle != null) {
    return (
      getHandleType(connection.sourceHandle) ===
      getHandleType(connection.targetHandle)
    )
  } else {
    return true
  }
}

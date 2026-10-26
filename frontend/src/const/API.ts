const HOST =
  process.env.REACT_APP_SERVER_HOST || window.location.hostname || "localhost"
const PROTO =
  process.env.REACT_APP_SERVER_PROTO ||
  window.location.protocol.replace(":", "") ||
  "http"
const PORT =
  process.env.REACT_APP_SERVER_PORT ||
  window.location.port ||
  (PROTO == "https" ? 443 : 80)

export const BASE_URL =
  PORT == null ? `${PROTO}://${HOST}` : `${PROTO}://${HOST}:${PORT}`

/**
 * API Request Timeout Constants (in milliseconds, set to axios)
 */
export const API_TIMEOUT = {
  DEFAULT: 600000, // Default timeout for backend API
  UPLOAD_DOWNLOAD: 0, // Upload/Download API timeout (No timeout(0), consider large file transfers)
} as const

import type {
  ExperimentDTO,
  ExperimentsDTO,
  FunctionsDTO,
} from "api/experiments/Experiments"
import { RunResultDTO } from "api/run/Run"
import type {
  ExperimentListType,
  ExperimentType,
  ExperimentFunction,
  EXPERIMENTS_STATUS,
} from "store/slice/Experiments/ExperimentsType"

export function convertToExperimentListType(
  dto: ExperimentsDTO,
): ExperimentListType {
  const experimentList: ExperimentListType = {}
  Object.entries(dto).forEach(([uid, value]) => {
    experimentList[uid] = convertToExperimentType(value)
  })
  return experimentList
}

export function convertToExperimentType(dto: ExperimentDTO): ExperimentType {
  const functions: { [nodeId: string]: ExperimentFunction } = {}
  Object.entries(dto.function).forEach(([_, value]) => {
    const status = convertToExperimentStatus(value.success)
    functions[value.unique_id] = {
      name: value.name,
      nodeId: value.unique_id,
      status: status,
      hasNWB: value.hasNWB,
    }
    status === "error" &&
      value.message &&
      (functions[value.unique_id].message = value.message)
  })
  return {
    uid: dto.unique_id,
    startedAt: dto.started_at,
    finishedAt: dto.finished_at,
    status: dto.success,
    name: dto.name,
    hasNWB: dto.hasNWB,
    functions,
    frameRate: dto.nwb?.imaging_plane.imaging_rate,
  }
}

function convertToExperimentStatus(dto: string): EXPERIMENTS_STATUS {
  switch (dto) {
    case "running":
      return "running"
    case "success":
      return "success"
    case "error":
      return "error"
    default:
      throw new Error("failed to convert to EXPERIMENTS_STATUS")
  }
}

export function convertFunctionsToRunResultDTO(
  dto: FunctionsDTO,
): RunResultDTO {
  const result: RunResultDTO = {}
  Object.entries(dto).forEach(([nodeId, value]) => {
    result[nodeId] = {
      status: value.success,
      message: value.message ?? "",
      name: value.name,
      outputPaths: value.outputPaths,
    }
  })
  return result
}

export function convertToFlowChartList(
  dto: ExperimentsDTO,
): ExperimentListType {
  const experimentList: ExperimentListType = {}
  Object.entries(dto).forEach(([uid, value]) => {
    experimentList[uid] = convertToExperimentType(value)
  })
  return experimentList
}

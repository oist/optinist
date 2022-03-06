import type { ExperimentDTO, ExperimentsDTO } from 'api/experiments/Experiments'
import type {
  ExperimentListType,
  ExperimentType,
  ExperimentFunction,
  EXPERIMENTS_STATUS,
} from './ExperimentsType'

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
  Object.entries(dto.function).forEach(([name, value]) => {
    functions[value.unique_id] = {
      name,
      nodeId: value.unique_id,
      postion: value.position,
      status: convertToExperimentStatus(value.success),
    }
  })
  return {
    uid: dto.unique_id,
    timestamp: dto.timestamp,
    name: dto.name,
    functions,
  }
}

function convertToExperimentStatus(dto: string): EXPERIMENTS_STATUS {
  switch (dto) {
    case 'running':
      return 'running'
    case 'success':
      return 'success'
    case 'error':
      return 'error'
    default:
      throw new Error('failed to convert to EXPERIMENTS_STATUS')
  }
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

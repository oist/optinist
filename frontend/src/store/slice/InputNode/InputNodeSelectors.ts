import { getAllFileTypeConfigs, FILE_TYPE } from "config/fileTypes.config"
import { FileNodeFactory } from "factories/FileNodeFactory"
import { InputNodeType } from "store/slice/InputNode/InputNodeType"
import {
  isImageInputNode,
  isCsvInputNode,
  isHDF5InputNode,
  isMatlabInputNode,
  isMicroscopeInputNode,
} from "store/slice/InputNode/InputNodeUtils"
import { RootState } from "store/store"

// Generic selector factory for typed file path selectors
function createTypedFilePathSelector(
  fileType: FILE_TYPE,
  typePredicate: (node: InputNodeType) => boolean,
) {
  return (nodeId: string) => (state: RootState) => {
    const node = selectInputNodeById(nodeId)(state)
    if (typePredicate(node)) {
      const config = FileNodeFactory.getFileTypeConfig(fileType)
      if (config?.filePathType === "array") {
        return node.selectedFilePath as string[] | undefined
      } else {
        return node.selectedFilePath as string | undefined
      }
    } else {
      throw new Error("invalid input node type")
    }
  }
}

export const selectInputNode = (state: RootState) => state.inputNode

export const selectInputNodeById = (nodeId: string) => (state: RootState) =>
  state.inputNode[nodeId]

export const selectInputNodeDefined = (nodeId: string) => (state: RootState) =>
  Object.keys(state.inputNode).includes(nodeId)

export const selectInputNodeFileType = (nodeId: string) => (state: RootState) =>
  selectInputNodeById(nodeId)(state).fileType

export const selectInputNodeSelectedFilePath =
  (nodeId: string) => (state: RootState) => {
    return selectInputNodeById(nodeId)(state).selectedFilePath
  }

// Dynamically generate typed selectors from config
const generateTypedSelectors = () => {
  const selectors: Record<
    string,
    (nodeId: string) => (state: RootState) => string | string[] | undefined
  > = {}

  // Map file types to their predicates
  const typePredicateMap: Record<FILE_TYPE, (node: InputNodeType) => boolean> =
    {
      image: isImageInputNode,
      csv: isCsvInputNode,
      hdf5: isHDF5InputNode,
      matlab: isMatlabInputNode,
      microscope: isMicroscopeInputNode,
    }

  getAllFileTypeConfigs().forEach((config) => {
    const selectorName = FileNodeFactory.generateSelectorName(
      config.key,
      "SelectedFilePath",
    )
    const predicate = typePredicateMap[config.key]
    if (predicate) {
      selectors[selectorName] = createTypedFilePathSelector(
        config.key,
        predicate,
      )
    }
  })

  return selectors
}

const typedSelectors = generateTypedSelectors()

// Export dynamically generated selectors
export const selectCsvInputNodeSelectedFilePath =
  typedSelectors.selectCsvInputNodeSelectedFilePath
export const selectImageInputNodeSelectedFilePath =
  typedSelectors.selectImageInputNodeSelectedFilePath
export const selectHdf5InputNodeSelectedFilePath =
  typedSelectors.selectHdf5InputNodeSelectedFilePath
export const selectHDF5InputNodeSelectedFilePath =
  typedSelectors.selectHdf5InputNodeSelectedFilePath // Backward compatibility
export const selectMatlabInputNodeSelectedFilePath =
  typedSelectors.selectMatlabInputNodeSelectedFilePath
export const selectMicroscopeInputNodeSelectedFilePath =
  typedSelectors.selectMicroscopeInputNodeSelectedFilePath

// Dynamic selectors generation capability is available for future extensions via getAllFileTypeConfigs()

export const selectFilePathIsUndefined = (state: RootState) =>
  Object.keys(state.inputNode).length === 0 ||
  Object.values(state.inputNode).filter((inputNode) => {
    if (isHDF5InputNode(inputNode)) {
      return inputNode.selectedFilePath == null || inputNode.hdf5Path == null
    } else if (isMatlabInputNode(inputNode)) {
      return inputNode.selectedFilePath == null || inputNode.matPath == null
    } else {
      const filePath = inputNode.selectedFilePath
      if (Array.isArray(filePath)) {
        return filePath.length === 0
      } else {
        return filePath == null
      }
    }
  }).length > 0

export const selectInputNodeParam = (nodeId: string) => (state: RootState) =>
  selectInputNodeById(nodeId)(state).param

const selectCsvInputNodeParam = (nodeId: string) => (state: RootState) => {
  const inputNode = selectInputNodeById(nodeId)(state)
  if (isCsvInputNode(inputNode)) {
    return inputNode.param
  } else {
    throw new Error(`The InputNode is not CsvInputNode. (nodeId: ${nodeId})`)
  }
}

export const selectCsvInputNodeParamSetHeader =
  (nodeId: string) => (state: RootState) =>
    selectCsvInputNodeParam(nodeId)(state).setHeader

export const selectCsvInputNodeParamSetIndex =
  (nodeId: string) => (state: RootState) =>
    selectCsvInputNodeParam(nodeId)(state).setIndex

export const selectCsvInputNodeParamTranspose =
  (nodeId: string) => (state: RootState) =>
    selectCsvInputNodeParam(nodeId)(state).transpose

export const selectInputNodeHDF5Path =
  (nodeId: string) => (state: RootState) => {
    const item = selectInputNodeById(nodeId)(state)
    if (isHDF5InputNode(item)) {
      return item.hdf5Path
    } else {
      return undefined
    }
  }

export const selectInputNodeMatlabPath =
  (nodeId: string) => (state: RootState) => {
    const item = selectInputNodeById(nodeId)(state)
    if (isMatlabInputNode(item)) {
      return item.matPath
    } else {
      return undefined
    }
  }

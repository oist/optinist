import {
  configureStore,
  ThunkAction,
  Action,
  combineReducers,
} from '@reduxjs/toolkit'
import {Elements} from 'react-flow-renderer'
import {
  algorithmListReducer,
  algorithmNodeReducer,
  displayDataReducer,
  fileUploaderReducer,
  flowElementReducer,
  inputNodeReducer,
  nwbReducer,
  filesTreeReducer,
  handleTypeColorReducer,
  rightDrawerReducer,
  visualaizeItemReducer,
  snakemakeReducer,
  pipelineReducer,
  hdf5Reducer,
  experimentsReducer,
} from './slice'
import {AlgorithmListType} from "./slice/AlgorithmList/AlgorithmListType";
import {AlgorithmNode} from "./slice/AlgorithmNode/AlgorithmNodeType";
import {DisplayData} from "./slice/DisplayData/DisplayDataType";
import {FileUploader} from "./slice/FileUploader/FileUploaderType";
import {NodeData} from "./slice/FlowElement/FlowElementType";
import {InputNode} from "./slice/InputNode/InputNodeType";
import {HandleTypeColor} from "./slice/HandleTypeColor/HandleTypeColorSlice";
import {FilesTree} from "./slice/FilesTree/FilesTreeType";
import {NWBType} from "./slice/NWB/NWBType";
import {RIGHT_DRAWER_MODE_TYPE} from "./slice/RightDrawer/RightDrawerSlice";
import {VisualaizeItem} from "./slice/VisualizeItem/VisualizeItemType";
import {SnakemakeType} from "./slice/Snakemake/SnakemakeType";
import {Pipeline} from "./slice/Pipeline/PipelineType";
import {HDF5Tree} from "./slice/HDF5/HDF5Type";
import {Experiments} from "./slice/Experiments/ExperimentsType";

export type State = {
  algorithmList: AlgorithmListType
  algorithmNode: AlgorithmNode
  displayData: DisplayData
  fileUploader: FileUploader
  flowElement: Elements<NodeData>
  inputNode: InputNode
  handleColor: HandleTypeColor
  filesTree: FilesTree
  nwb: NWBType
  rightDrawer: RIGHT_DRAWER_MODE_TYPE
  visualaizeItem: VisualaizeItem
  snakemake: SnakemakeType
  pipeline: Pipeline
  hdf5: HDF5Tree
  experiments: Experiments
}
export const rootReducer = combineReducers({
  algorithmList: algorithmListReducer,
  algorithmNode: algorithmNodeReducer,
  displayData: displayDataReducer,
  fileUploader: fileUploaderReducer,
  flowElement: flowElementReducer,
  inputNode: inputNodeReducer,
  handleColor: handleTypeColorReducer,
  filesTree: filesTreeReducer,
  nwb: nwbReducer,
  rightDrawer: rightDrawerReducer,
  visualaizeItem: visualaizeItemReducer,
  snakemake: snakemakeReducer,
  pipeline: pipelineReducer,
  hdf5: hdf5Reducer,
  experiments: experimentsReducer,
})

export const store = configureStore({
  reducer: rootReducer,
})

export type AppDispatch = typeof store.dispatch
export type RootState = ReturnType<typeof store.getState>
export type AppThunk<ReturnType = void> = ThunkAction<
  ReturnType,
  RootState,
  unknown,
  Action<string>
>
export type ThunkApiConfig<T = unknown, PendingMeta = unknown> = {
  state: RootState
  dispatch: AppDispatch
  rejectValue: T
  penfingMeta: PendingMeta
}

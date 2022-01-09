import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { Layout, Model, TabNode } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import './App.css'
import { flexjson } from 'const/flexlayout'
// import { AlgorithmTreeView } from 'components/TreeView'
// import { FlowChart } from 'components/FlowChart'
import { ParamForm } from 'components/FlowChart/ParamForm'
import { Plot } from 'components/Plot'
import { ToolBar } from 'components/ToolBar'
import {
  selectDataTypeByTabId,
  selectFilePathByTabId,
  selectNodeIdByTabId,
  selectTabIsList,
} from 'store/slice/LayoutTab/LayoutTabSelectors'
import { arrayEqualityFn } from 'utils/EqualityUtils'
import { deleteAllTabsByIdList } from 'store/slice/LayoutTab/LayputTabSlice'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import { TAB_COMPONENT_TYPE_SET } from 'store/slice/LayoutTab/LayoutTabType'
import AppLayout from './components/Layout'
const model = Model.fromJson(flexjson)

export const FlexLayoutModelContext = React.createContext<Model>(model)
/**
 * nodeId
 */
export const ParamFormContext = React.createContext<string>('') // todo 後で場所移動
export const DisplayDataTabContext = React.createContext<{
  nodeId: string
  filePath: string
  dataType: DATA_TYPE
}>({ nodeId: '', filePath: '', dataType: 'table' })

// const Factory = (node: TabNode) => {
//   var component = node.getComponent()
//   const layoutTabId = node.getId()
//   switch (component) {
//     case 'flowchart':
//       return <FlowChart />
//     case 'sidebar':
//       return <SideBar />
//     case TAB_COMPONENT_TYPE_SET.PARAM_FORM:
//       return <ParamFormTab layoutTabId={layoutTabId} />
//     case TAB_COMPONENT_TYPE_SET.DISPLAY_DATA:
//       return <DisplayTab layoutTabId={layoutTabId} />
//     default:
//       return null
//   }
// }

// function App() {
//   // modelのtabが削除された場合にstateの方も削除する
//   const dispatch = useDispatch()
//   const stateTabIdList = useSelector(selectTabIsList, arrayEqualityFn)
//   React.useEffect(() => {
//     const remainderTabIdList = stateTabIdList.filter(
//       (tabId) => model.getNodeById(tabId) === null,
//     )
//     if (remainderTabIdList.length > 0) {
//       dispatch(deleteAllTabsByIdList(remainderTabIdList))
//     }
//   }, [stateTabIdList, dispatch])

//   return (
//     <div id="container">
//       <div className="app">
//         <FlexLayoutModelContext.Provider value={model}>
//           <div className="toolbar">
//             <ToolBar />
//           </div>
//           <div className="contents">
//             <Layout model={model} factory={Factory} />
//           </div>
//         </FlexLayoutModelContext.Provider>
//       </div>
//     </div>
//   )
// }

// const DisplayTab = React.memo<{ layoutTabId: string }>(({ layoutTabId }) => {
//   const nodeId = useSelector(selectNodeIdByTabId(layoutTabId))
//   const dataType = useSelector(selectDataTypeByTabId(layoutTabId))
//   const filePath = useSelector(selectFilePathByTabId(layoutTabId))
//   if (filePath != null && dataType != null) {
//     return (
//       <DisplayDataTabContext.Provider value={{ nodeId, filePath, dataType }}>
//         <Plot />
//       </DisplayDataTabContext.Provider>
//     )
//   } else {
//     return null
//   }
// })

// const ParamFormTab = React.memo<{ layoutTabId: string }>(({ layoutTabId }) => {
//   const nodeId = useSelector(selectNodeIdByTabId(layoutTabId))
//   return (
//     <ParamFormTabContext.Provider value={nodeId}>
//       <ParamForm />
//     </ParamFormTabContext.Provider>
//   )
// })

const App: React.FC = () => {
  return <AppLayout />
}

export default App

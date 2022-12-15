import { store, rootReducer } from 'store/store'
import {
  pushInitialItemToNewRow,
  insertInitialItemToNextColumn,
  addItemForWorkflowDialog,
  deleteAllItemForWorkflowDialog,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  deleteDisplayItem,
  setNewDisplayDataPath,
} from 'store/slice/VisualizeItem/VisualizeItemActions'
import {
  selectVisualizeItemById,
  selectVisualizeItemLayout,
  selectSelectedVisualizeItemId,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { DATA_TYPE_SET } from 'store/slice/DisplayData/DisplayDataType'
import {
  isCsvItem,
  isImageItem,
} from 'store/slice/VisualizeItem/VisualizeItemUtils'
import { selectImageData } from 'store/slice/DisplayData/DisplayDataSelectors'

describe('Visualize', () => {
  const initialRootState = store.getState()

  // Item追加(行)
  test(pushInitialItemToNewRow.type, () => {
    // 1回目
    const targetState = rootReducer(initialRootState, pushInitialItemToNewRow)
    const itemId = 0
    expect(selectVisualizeItemById(itemId)(targetState)).toBeDefined()
    expect(selectVisualizeItemLayout(targetState)).toEqual([[itemId]])
    expect(selectSelectedVisualizeItemId(targetState)).toBe(itemId)
    expect(selectVisualizeItemById(itemId)(targetState).isWorkflowDialog).toBe(
      false,
    )
    // 2回目
    const nextTargetState = rootReducer(targetState, pushInitialItemToNewRow)
    const nextItemId = itemId + 1
    expect(selectVisualizeItemById(nextItemId)(nextTargetState)).toBeDefined()
    expect(selectVisualizeItemLayout(nextTargetState)).toEqual([
      [itemId],
      [nextItemId],
    ])
    expect(selectSelectedVisualizeItemId(nextTargetState)).toBe(nextItemId)
    expect(
      selectVisualizeItemById(itemId)(nextTargetState).isWorkflowDialog,
    ).toBe(false)
  })

  // Item追加(列)
  test(insertInitialItemToNextColumn.type, () => {
    const prevState = rootReducer(initialRootState, pushInitialItemToNewRow)
    const prevItemId = 0
    const targetState = rootReducer(
      prevState,
      insertInitialItemToNextColumn(prevItemId),
    )
    const itemId = prevItemId + 1
    expect(selectVisualizeItemById(itemId)(targetState)).toBeDefined()
    expect(selectVisualizeItemLayout(targetState)).toEqual([
      [prevItemId, itemId],
    ])
    expect(selectSelectedVisualizeItemId(targetState)).toBe(itemId)
    expect(selectVisualizeItemById(itemId)(targetState).isWorkflowDialog).toBe(
      false,
    )
  })

  // filePath変更

  test(`${setNewDisplayDataPath.type} when deleteData is false`, () => {
    const prevState = rootReducer(initialRootState, pushInitialItemToNewRow)
    const itemId = 0
    const nodeId = 'dummy_image2image8time_qikqvda3um'
    const filePath =
      '/tmp/studio/output/4fa05cff/dummy_image2image8time_qikqvda3um/image.json'
    const dataType = DATA_TYPE_SET.IMAGE
    const deleteData = false
    const targetState = rootReducer(
      prevState,
      setNewDisplayDataPath({
        itemId,
        nodeId,
        filePath,
        dataType,
        deleteData,
      }),
    )
    expect(isImageItem(selectVisualizeItemById(itemId)(targetState))).toBe(true)
    expect(selectVisualizeItemById(itemId)(targetState).nodeId).toBe(nodeId)
    expect(selectVisualizeItemById(itemId)(targetState).filePath).toBe(filePath)
  })

  test(`${setNewDisplayDataPath.type} when deleteData is true`, () => {
    const itemId = 0
    const prevNodeId = 'input_0'
    const prevFilePath = '/tmp/studio/input/hoge/hoge.tif'
    const prevDataType = DATA_TYPE_SET.IMAGE
    let prevState = rootReducer(
      rootReducer(initialRootState, pushInitialItemToNewRow),
      setNewDisplayDataPath({
        itemId,
        nodeId: prevNodeId,
        filePath: prevFilePath,
        dataType: prevDataType,
        deleteData: false,
      }),
    )
    prevState = {
      ...prevState,
      displayData: {
        ...prevState.displayData,
        image: {
          [prevFilePath]: {
            type: prevDataType,
            data: [[]],
            pending: false,
            fulfilled: true,
            error: null,
          },
        },
      },
    }
    const deleteData = true
    const filePath = '/tmp/studio/input/sample/sample.csv'
    const nodeId = 'input_90i2lotktx'
    const dataType = DATA_TYPE_SET.CSV
    const targetState = rootReducer(
      prevState,
      setNewDisplayDataPath({
        itemId,
        nodeId,
        filePath,
        dataType,
        deleteData,
        prevDataType,
        prevFilePath,
      }),
    )
    expect(isCsvItem(selectVisualizeItemById(itemId)(targetState))).toBe(true)
    expect(selectVisualizeItemById(itemId)(targetState).nodeId).toBe(nodeId)
    expect(selectVisualizeItemById(itemId)(targetState).filePath).toBe(filePath)
    // displayData削除
    expect(selectImageData(prevFilePath)(targetState)).toBeUndefined()
  })

  // Item削除

  test(`${deleteDisplayItem.type} when deleteData is false`, () => {
    const prevState = rootReducer(initialRootState, pushInitialItemToNewRow)
    const itemId = 0
    const targetState = rootReducer(
      prevState,
      deleteDisplayItem({ itemId, deleteData: false }),
    )
    expect(selectVisualizeItemById(itemId)(targetState)).toBeUndefined()
    expect(selectVisualizeItemLayout(targetState).flat().includes(itemId)).toBe(
      false,
    )
  })

  test(`${deleteDisplayItem.type} when deleteData is true`, () => {
    const itemId = 0
    const prevNodeId = 'input_0'
    const prevFilePath = '/tmp/studio/input/hoge/hoge.tif'
    const prevDataType = DATA_TYPE_SET.IMAGE
    let prevState = rootReducer(
      rootReducer(initialRootState, pushInitialItemToNewRow),
      setNewDisplayDataPath({
        itemId,
        nodeId: prevNodeId,
        filePath: prevFilePath,
        dataType: prevDataType,
        deleteData: false,
      }),
    )
    prevState = {
      ...prevState,
      displayData: {
        ...prevState.displayData,
        image: {
          [prevFilePath]: {
            type: prevDataType,
            data: [[]],
            pending: false,
            fulfilled: true,
            error: null,
          },
        },
      },
    }
    const deleteData = true
    const targetState = rootReducer(
      prevState,
      deleteDisplayItem({
        itemId,
        deleteData,
        filePath: prevFilePath,
        dataType: prevDataType,
      }),
    )
    expect(selectVisualizeItemById(itemId)(targetState)).toBeUndefined()
    expect(selectVisualizeItemLayout(targetState).flat().includes(itemId)).toBe(
      false,
    )
    // displayData削除
    expect(selectImageData(prevFilePath)(targetState)).toBeUndefined()
  })

  // AlgorithmNodeのoutputダイアログで使用するItemの追加
  test(addItemForWorkflowDialog.type, () => {
    const itemId = 0
    const nodeId = 'dummy_image2image8time_qikqvda3um'
    const filePath =
      '/tmp/studio/output/4fa05cff/dummy_image2image8time_qikqvda3um/image.json'
    const dataType = DATA_TYPE_SET.IMAGE
    const targetState = rootReducer(
      initialRootState,
      addItemForWorkflowDialog({ nodeId, filePath, dataType }),
    )
    expect(isImageItem(selectVisualizeItemById(itemId)(targetState))).toBe(true)
    expect(selectVisualizeItemById(itemId)(targetState).isWorkflowDialog).toBe(
      true,
    )
    expect(selectVisualizeItemById(itemId)(targetState).nodeId).toBe(nodeId)
    expect(selectVisualizeItemById(itemId)(targetState).filePath).toBe(filePath)
    expect(selectVisualizeItemLayout(targetState).flat().includes(itemId)).toBe(
      false,
    )
  })

  // AlgorithmNodeのoutputダイアログを閉じる
  test(deleteAllItemForWorkflowDialog.type, () => {
    const nodeId = 'dummy_image2image8time_qikqvda3um'
    // itemId: 0
    let prevState = rootReducer(
      initialRootState,
      addItemForWorkflowDialog({
        nodeId,
        filePath:
          '/tmp/studio/output/4fa05cff/dummy_image2image8time_qikqvda3um/image.json',
        dataType: DATA_TYPE_SET.IMAGE,
      }),
    )
    // itemId: 1
    prevState = rootReducer(
      prevState,
      addItemForWorkflowDialog({
        nodeId,
        filePath:
          '/tmp/studio/output/4fa05cff/dummy_image2image8time_qikqvda3um/timeseries',
        dataType: DATA_TYPE_SET.TIME_SERIES,
      }),
    )
    const targetState = rootReducer(prevState, deleteAllItemForWorkflowDialog())
    expect(selectVisualizeItemById(0)(targetState)).toBeUndefined()
    expect(selectVisualizeItemById(1)(targetState)).toBeUndefined()
  })
})

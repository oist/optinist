import { TreeNodeTypeDTO } from "api/files/Files"
import { getFilesTree } from "store/slice/FilesTree/FilesTreeAction"
import reducer, { initialState } from "store/slice/FilesTree/FilesTreeSlice"
import { FilesTree } from "store/slice/FilesTree/FilesTreeType"

describe("FilesTree", () => {
  const mockPayload: TreeNodeTypeDTO[] = [
    {
      path: "/tmp/optinist/input/hoge",
      name: "hoge",
      isdir: true,
      nodes: [
        {
          path: "/tmp/optinist/input/hoge/hoge.tif",
          name: "hoge.tif",
          isdir: false,
        },
      ],
    },
    {
      path: "/tmp/optinist/input/copy_image1",
      name: "copy_image1",
      isdir: true,
      nodes: [
        {
          path: "/tmp/optinist/input/copy_image1/copy_image1.tif",
          name: "copy_image1.tif",
          isdir: false,
        },
      ],
    },
  ]

  const expectState: FilesTree = {
    image: {
      isLoading: false,
      isLatest: true,
      tree: [
        {
          path: "/tmp/optinist/input/hoge",
          name: "hoge",
          isDir: true,
          nodes: [
            {
              path: "/tmp/optinist/input/hoge/hoge.tif",
              name: "hoge.tif",
              isDir: false,
            },
          ],
        },
        {
          path: "/tmp/optinist/input/copy_image1",
          name: "copy_image1",
          isDir: true,
          nodes: [
            {
              path: "/tmp/optinist/input/copy_image1/copy_image1.tif",
              name: "copy_image1.tif",
              isDir: false,
            },
          ],
        },
      ],
    },
  }

  test(getFilesTree.fulfilled.type, () => {
    expect(
      reducer(
        reducer(initialState, {
          type: getFilesTree.pending.type,
          meta: {
            arg: { fileType: "image" },
            requestId: "F0QeIMS-KV132B2q79qaz",
            requestStatus: "pending",
          },
        }),
        {
          type: getFilesTree.fulfilled.type,
          payload: mockPayload,
          meta: {
            arg: { fileType: "image" },
            requestId: "F0QeIMS-KV132B2q79qaz",
            requestStatus: "fulfilled",
          },
        },
      ),
    ).toEqual(expectState)
  })
})

export const mockStoreData = {
  isLatest: true,
  algorithmList: {
    tree: {
      caiman: {
        type: "parent",
        children: {
          caiman_mc: {
            type: "child",
            functionPath: "caiman/caiman_mc",
            args: [{ name: "image", type: "ImageData", isNone: false }],
            returns: [{ name: "mc_images", type: "ImageData" }],
          },
          caiman_cnmf: {
            type: "child",
            functionPath: "caiman/caiman_cnmf",
            args: [{ name: "images", type: "ImageData", isNone: false }],
            returns: [
              { name: "fluorescence", type: "FluoData" },
              { name: "iscell", type: "IscellData" },
            ],
          },
        },
      },
      suite2p: {
        type: "parent",
        children: {
          suite2p_file_convert: {
            type: "child",
            functionPath: "suite2p/suite2p_file_convert",
            args: [{ name: "image", type: "ImageData", isNone: false }],
            returns: [{ name: "ops", type: "Suite2pData" }],
          },
          suite2p_registration: {
            type: "child",
            functionPath: "suite2p/suite2p_registration",
            args: [{ name: "ops", type: "Suite2pData", isNone: false }],
            returns: [{ name: "ops", type: "Suite2pData" }],
          },
          // 残りのデータも同様に記述
        },
      },
      lccd: {
        type: "parent",
        children: {
          lccd_file_convert: {
            type: "child",
            functionPath: "lccd/lccd_file_convert",
            args: [{ name: "image", type: "ImageData", isNone: false }],
            returns: [{ name: "lccd", type: "LccdData" }],
          },
          // 残りのデータも同様に記述
        },
      },
      optinist: {
        type: "parent",
        children: {
          optinist_file_convert: {
            type: "child",
            functionPath: "optinist/optinist_file_convert",
            args: [{ name: "image", type: "ImageData", isNone: false }],
            returns: [{ name: "optinist", type: "OptinistData" }],
          },
          // 残りのデータも同様に記述
        },
      },
    },
  },
  pipeline: {
    currentPipeline: {
      uid: "123",
      name: "pipeline1",
    },
  },
}

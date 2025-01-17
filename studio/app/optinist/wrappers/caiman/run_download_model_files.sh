#!/bin/bash

#
# This function is a port of
#   studio.app.optinist.wrappers.caiman.cnmf.util_download_model_files()
#
download_model_files() {
  # Base URL for downloading model files
  BASE_URL="https://github.com/flatironinstitute/CaImAn/raw/v1.9.12/model"
  
  # List of model files to download
  MODEL_FILES=(
      "cnn_model.h5"
      "cnn_model.h5.pb"
      "cnn_model.json"
      "cnn_model_online.h5"
      "cnn_model_online.h5.pb"
      "cnn_model_online.json"
  )
  
  # Create caiman_data directory in home directory
  CAIMAN_DATA_DIR="$HOME/caiman_data"
  if [ ! -d "$CAIMAN_DATA_DIR" ]; then
    mkdir -p "$CAIMAN_DATA_DIR"
  fi
  
  # Create model directory
  MODEL_DIR="$CAIMAN_DATA_DIR/model"
  if [ ! -d "$MODEL_DIR" ]; then
    mkdir -p "$MODEL_DIR"
  fi
  
  # Count existing files in model directory
  FILE_COUNT=$(ls -1 "$MODEL_DIR" | wc -l)
  
  # Download files if needed
  if [ "$FILE_COUNT" -lt "${#MODEL_FILES[@]}" ]; then
      for MODEL in "${MODEL_FILES[@]}"; do
          FILE_PATH="$MODEL_DIR/$MODEL"
          if [ ! -f "$FILE_PATH" ]; then
              echo "Downloading $MODEL"
              curl -L "$BASE_URL/$MODEL" -o "$FILE_PATH"
          fi
      done
  fi
}

# call downloading func
download_model_files

import os

# set environment variables for development.
os.environ["OPTINIST_DEV_ROOT_DIR"] = os.path.dirname(os.path.abspath(__file__))

# run backend main module.
if __name__ == "__main__":
    from optinist.__main_unit__ import main
    main(develop_mode = True)

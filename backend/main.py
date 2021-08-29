from fastapi import Depends, FastAPI
import uvicorn

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.get("/params")
async def params():
    import params
    return params.get_params()

@app.get("/run")
async def run():
    import run
    return run.run_code()


if __name__ == '__main__':
	uvicorn.run(app, host="0.0.0.0", port=8000)

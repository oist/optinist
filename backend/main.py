from fastapi import Depends, FastAPI
import uvicorn
from starlette.middleware.cors import CORSMiddleware

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.get("/params/{name}")
async def params(name: str):
    if name == 'test':
        import params
        return params.get_params()

@app.get("/run")
async def run():
    import run
    return run.run_code()


if __name__ == '__main__':
	uvicorn.run(app, host="0.0.0.0", port=8000)

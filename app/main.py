from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.routes import router

app = FastAPI(
    title="BioTech Genomic Sequence Analytics Pipeline",
    description="Algorithmic bioinformatics engine leveraging Python modules for nucleotide sequence transcription, protein translation, variant SNP detection, and pairwise alignment scoring.",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router)

@app.get("/health", tags=["System"])
def health_check():
    return {
        "status": "HEALTHY",
        "service": "BioTech Genomic Sequence Analytics Pipeline",
        "pipeline_version": "1.0.0",
        "supported_codecs": ["DNA", "RNA", "PROTEIN"]
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="0.0.0.0", port=8002, reload=True)

import os
from dotenv import load_dotenv

# Load environment variables from .env file if it exists
load_dotenv()

class Config:
    LM_STUDIO_API_BASE = os.getenv("LM_STUDIO_API_BASE", "http://localhost:1234/v1")
    LM_STUDIO_API_KEY = os.getenv("LM_STUDIO_API_KEY", "lm-studio")
    ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "your.email@example.com")

if Config.ENTREZ_EMAIL == "your.email@example.com":
    print("WARNING: Please set your specific email in .env for PubMed/Entrez access.")

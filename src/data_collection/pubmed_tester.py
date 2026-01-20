import sys
import os

# Add the project root to sys.path to allow importing from src
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from Bio import Entrez
from src.config import Config

def test_pubmed_query(query):
    """
    Tests a PubMed query using Bio.Entrez and returns the count of results.
    """
    Entrez.email = Config.ENTREZ_EMAIL
    
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
        record = Entrez.read(handle)
        handle.close()
        
        count = int(record["Count"])
        return count, None
    except Exception as e:
        return 0, str(e)

if __name__ == "__main__":
    query = "machine learning[Title]"
    print(f"Testing PubMed query: '{query}'")
    count, error = test_pubmed_query(query)
    
    if error:
        print(f"Error: {error}")
    else:
        print(f"Success! Found {count} results.")

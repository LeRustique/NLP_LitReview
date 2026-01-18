import pandas as pd
from Bio import Entrez
from src.config import Config
import sys

def repair_missing_dois(input_file):
    """
    Reads the input CSV, identifies rows with missing DOIs but having a PMID/ID.
    Fetches the DOI from PubMed and updates the CSV.
    """
    print(f"--- Repairing DOIs for {input_file} ---")
    
    Entrez.email = Config.ENTREZ_EMAIL
    
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        return

    # Add DOI column if missing
    if "DOI" not in df.columns:
        df["DOI"] = ""

    # Find rows to update
    # Logic: Source=PubMed or ID_Type contains PMID, and DOI is empty
    mask = ( (df["Source"] == "PubMed") | (df["ID_Type"].astype(str).str.contains("PMID")) ) & \
           ( (df["DOI"].isnull()) | (df["DOI"] == "") | (df["DOI"] == "nan") )
    
    rows_to_check = df[mask]
    
    if rows_to_check.empty:
        print("No PubMed articles with missing DOIs found.")
        return

    print(f"Found {len(rows_to_check)} PubMed articles missing DOI. Fetching...")

    # We can batch fetch using PMIDs
    # Extract PMIDs (assuming 'ID' column holds the PMID)
    pmids = []
    # Map PMID -> (index in df)
    pmid_map = {}
    
    for index, row in rows_to_check.iterrows():
        pmid = str(row["ID"]).strip()
        if pmid.isdigit():
            pmids.append(pmid)
            if pmid not in pmid_map:
                pmid_map[pmid] = []
            pmid_map[pmid].append(index)
            
    if not pmids:
        print("No valid PMIDs found in the missing rows.")
        return

    # Fetch in batches of 200
    batch_size = 200
    updated_count = 0
    
    for i in range(0, len(pmids), batch_size):
        batch_pmids = pmids[i:i+batch_size]
        print(f"  Fetching batch {i+1}-{i+len(batch_pmids)}...")
        
        try:
            handle = Entrez.efetch(db="pubmed", id=",".join(batch_pmids), retmode="xml")
            data = Entrez.read(handle)
            handle.close()
            
            if 'PubmedArticle' in data:
                for article in data['PubmedArticle']:
                    medline = article.get("MedlineCitation", {})
                    pmid_res = str(medline.get("PMID", ""))
                    
                    # Extract DOI
                    doi = ""
                    article_data = medline.get("Article", {})
                    
                    # Try ELocationID
                    if "ELocationID" in article_data:
                        for eid in article_data["ELocationID"]:
                            if eid.attributes.get("EIdType") == "doi":
                                doi = str(eid)
                                break
                    
                    # Try PubmedData (preferred)
                    pubmed_data = article.get("PubmedData", {})
                    if "ArticleIdList" in pubmed_data:
                        for aid in pubmed_data["ArticleIdList"]:
                            if aid.attributes.get("IdType") == "doi":
                                doi = str(aid)
                                break
                    
                    if doi and pmid_res in pmid_map:
                        for idx in pmid_map[pmid_res]:
                            df.at[idx, "DOI"] = doi
                        updated_count += 1
                        
        except Exception as e:
            print(f"  Error fetching batch: {e}")

    print(f"Repair complete. Updated {updated_count} articles with DOIs.")
    
    # Save back
    try:
        df.to_csv(input_file, index=False, encoding='utf-8-sig')
        print(f"Saved updated file to {input_file}")
    except Exception as e:
        print(f"Error saving file: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        repair_missing_dois(sys.argv[1])
    else:
        # Default
        repair_missing_dois("results/final_selection.csv")

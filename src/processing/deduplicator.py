import pandas as pd
from thefuzz import fuzz
import re

def normalize_title(title):
    """
    Normalizes a title for comparison: lowercase, remove non-alphanumeric.
    """
    if not isinstance(title, str):
        return ""
    # Lowercase
    t = title.lower()
    # Remove punctuation/special chars
    t = re.sub(r'[^a-z0-9\s]', '', t)
    # Collapse whitespace
    t = re.sub(r'\s+', ' ', t).strip()
    return t

def deduplicate_data(pubmed_file, scholar_file, output_file):
    """
    Merges PubMed and Google Scholar data, removing duplicates.
    PubMed data is prioritized.
    """
    print("--- Starting Deduplication ---")
    
    # Load data
    try:
        df_pubmed = pd.read_csv(pubmed_file)
        print(f"Loaded {len(df_pubmed)} PubMed records.")
    except Exception:
        df_pubmed = pd.DataFrame()
        print("PubMed file not found or empty.")

    try:
        df_scholar = pd.read_csv(scholar_file)
        print(f"Loaded {len(df_scholar)} Scholar records.")
    except Exception:
        df_scholar = pd.DataFrame()
        print("Scholar file not found or empty.")

    if df_pubmed.empty and df_scholar.empty:
        print("No data to deduplicate.")
        return

    # Add Source if missing (for safety)
    if not df_pubmed.empty and "Source" not in df_pubmed.columns:
        df_pubmed["Source"] = "PubMed"
    if not df_scholar.empty and "Source" not in df_scholar.columns:
        df_scholar["Source"] = "Scholar"

    # Start with all PubMed records as unique
    unique_records = df_pubmed.to_dict('records')
    
    # Create a lookup for existing normalized titles
    # We compare against the list. 
    # For speed, we can use a set of normalized titles, but fuzzy match needs iteration.
    # Let's use exact normalized match for speed, and fuzzy for close calls if needed.
    
    existing_titles = set()
    for r in unique_records:
        t = normalize_title(r.get("Title", ""))
        if t:
            existing_titles.add(t)

    duplicates_found = 0
    new_records = 0
    
    scholar_records = df_scholar.to_dict('records')
    
    for rec in scholar_records:
        title = rec.get("Title", "")
        norm_title = normalize_title(title)
        
        if not norm_title:
             continue
             
        # 1. Exact Match on Normalized Title
        if norm_title in existing_titles:
            duplicates_found += 1
            # Optional: Merge info? E.g. Citations from Scholar -> PubMed record
            # Finding the original record is O(N) unless we built a dict.
            # For simplicity, we skip the Scholar duplicate.
            continue
            
        # 2. Fuzzy Match (Optional - can be slow O(N*M))
        # Only check if we are paranoid about slight differences.
        # Given "The" user request "eliminate duplicated", exact normalized is usually safe enough 
        # for "Perception de la Neuroréanimation" vs "Perception de la neuroréanimation".
        # Let's trust normalized match for now to be fast. 
        # If user wants fuzzy, we can enable it.
        
        # Add to unique
        unique_records.append(rec)
        existing_titles.add(norm_title)
        new_records += 1

    print(f"Deduplication complete.")
    print(f"  - PubMed entries kept: {len(df_pubmed)}")
    print(f"  - Scholar duplicates removed: {duplicates_found}")
    print(f"  - Scholar unique entries added: {new_records}")
    print(f"  - Total unique articles: {len(unique_records)}")

    # Save
    df_final = pd.DataFrame(unique_records)
    df_final.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"Saved merged dataset to {output_file}")

import os
from src.pubmed_tester import test_pubmed_query
from src.scholar_tester import test_scholar_query
from src.data_collection.pubmed_fetcher import fetch_pubmed_results
from src.data_collection.scholar_fetcher import fetch_scholar_results
from src.data_collection.data_manager import save_to_csv
from src.processing.deduplicator import deduplicate_data
from src.nlp_processor import NLPProcessor

def read_query_file(file_path):
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found.")
        return None
    with open(file_path, "r", encoding="utf-8") as f:
        return f.read().strip()

def main():
    print("--- Medical Literature Review Query Tester ---\n")
    
    # Paths to manual query files
    pubmed_path = os.path.join("request", "request_pubmed.md")
    scholar_path = os.path.join("request", "request_googleSchola.md")

    # 1. Test PubMed
    print(f"Reading PubMed query from {pubmed_path}...")
    pubmed_query = read_query_file(pubmed_path)
    
    if pubmed_query:
        print(f"Query: {pubmed_query}")
        print("Testing PubMed Query...")
        count, error = test_pubmed_query(pubmed_query)
        if error:
            print(f"  -> Error: {error}")
        else:
            print(f"  -> Found {count} results.")
    else:
        print("No PubMed query found.")

    print("")

    # 2. Test Google Scholar
    print(f"Reading Google Scholar query from {scholar_path}...")
    scholar_query = read_query_file(scholar_path)

    if scholar_query:
        print(f"Query: {scholar_query}")
        print("Testing Google Scholar Query...")
        count, val_text = test_scholar_query(scholar_query)
        if val_text and "Error" in val_text:
             print(f"  -> Error: {val_text}")
        elif val_text:
             print(f"  -> Found ~{count} results. ({val_text})")
        else:
             print(f"  -> Verification failed (could not parse page).")
    else:
        print("No Google Scholar query found.")

    print("\n--- Starting Data Collection ---\n")
    
    # 3. Collect PubMed Data
    if pubmed_query:
        print("Collecting data from PubMed...")
        pubmed_data = fetch_pubmed_results(pubmed_query)
        if pubmed_data:
            save_to_csv(pubmed_data, "results/pubmed_results.csv")
        else:
            print("No articles found to save.")
    
    # 4. Collect Google Scholar Data
    if scholar_query:
        print("\nCollecting data from Google Scholar (Max 200)...")
        print("Note: This performs scraping and may be slow or blocked.")
        scholar_data = fetch_scholar_results(scholar_query, max_results=200)
        if scholar_data:
            save_to_csv(scholar_data, "results/scholar_results.csv")
        else:
            print("No articles found or request blocked.")
            
    # 5. Deduplicate
    print("\n--- Processing Data ---\n")
    deduplicate_data(
        "results/pubmed_results.csv",
        "results/scholar_results.csv",
        "results/all_articles.csv"
    )

    # 6. NLP Screening
    print("\n--- Starting NLP Screening ---")
    
    # Load Criteria
    if not os.path.exists("criteria.md"):
        print("Warning: criteria.md not found. Skipping screening.")
        return
        
    with open("criteria.md", "r", encoding="utf-8") as f:
        criteria_text = f.read()

    # Load Data
    import pandas as pd
    try:
        df_articles = pd.read_csv("results/all_articles.csv")
    except Exception:
        print("Could not load results/all_articles.csv")
        return

    print(f"Screening {len(df_articles)} articles against criteria...")
    print("(This may take time depending on LLM speed)")

    processor = NLPProcessor()
    
    # Add columns if not exist
    if "Included" not in df_articles.columns:
        df_articles["Included"] = False
    if "Screening_Reason" not in df_articles.columns:
        df_articles["Screening_Reason"] = ""

    # Iterate (maybe limit to first 5 for testing?)
    # For production, iterate all.
    count = 0
    for index, row in df_articles.iterrows():
        # Check if already screened to allow resuming? 
        # For now, just screen all.
        
        # Skip if Abstract is empty/short
        if pd.isna(row["Abstract"]) or len(str(row["Abstract"])) < 20:
             df_articles.at[index, "Included"] = False
             df_articles.at[index, "Screening_Reason"] = "No Abstract"
             continue

        count += 1
        print(f"  Screening [{count}/{len(df_articles)}]: {str(row['Title'])[:50]}...")
        
        included, reason = processor.screen_paper(
            row["Title"], 
            row["Abstract"], 
            criteria_text
        )
        
        df_articles.at[index, "Included"] = included
        df_articles.at[index, "Screening_Reason"] = reason
        
        # Save periodically?
        if count % 10 == 0:
             df_articles.to_csv("results/screened_articles.csv", index=False, encoding='utf-8-sig')

    # Final Save
    df_articles.to_csv("results/screened_articles.csv", index=False, encoding='utf-8-sig')
    print(f"Screening complete. Results saved to results/screened_articles.csv")
    print(f"Included: {len(df_articles[df_articles['Included']==True])} papers.")

if __name__ == "__main__":
    main()

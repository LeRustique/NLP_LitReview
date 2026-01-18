import pandas as pd
import os
import sys

def manual_screen_articles(screened_file, output_file):
    """
    Interactively screens articles that were marked as 'Included' by the LLM.
    Updates the DataFrame with a 'Manual_Decision' column and saves it.
    """
    print("\n--- Starting Manual Screening ---\n")
    
    if not os.path.exists(screened_file):
        print(f"Error: {screened_file} not found.")
        return

    try:
        df = pd.read_csv(screened_file)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    if "Included" not in df.columns:
        print("Error: 'Included' column missing (LLM screening likely incomplete).")
        return

    # Initialize Manual_Decision column if not exists
    if "Manual_Decision" not in df.columns:
        df["Manual_Decision"] = "" # Empty strings by default

    # Columns to display
    display_cols = ["Title", "Abstract", "Journal", "Year", "Screening_Reason"]

    # Filter for candidates: Included by LLM AND not yet manually decided
    # We want to allow resuming, so we check if Manual_Decision is empty/null
    
    # Iterate through indices to update in-place
    total_candidates = len(df[df["Included"] == True])
    processed_count = 0
    
    for index, row in df.iterrows():
        # Only screen what LLM included
        if str(row["Included"]).lower() != "true":
            continue
            
        processed_count += 1
        
        # Skip if already decided (Resume capability)
        current_decision = str(row["Manual_Decision"]).lower()
        if current_decision in ["y", "n", "yes", "no", "included", "excluded"]:
            continue

        print("="*80)
        print(f"Paper {processed_count}/{total_candidates}")
        print("="*80)
        print(f"TITLE:    {row.get('Title', 'N/A')}")
        print("-" * 20)
        print(f"JOURNAL:  {row.get('Journal', 'N/A')} ({row.get('Year', 'N/A')})")
        print("-" * 20)
        print(f"LLM REASON: {row.get('Screening_Reason', 'N/A')}")
        print("-" * 20)
        print("ABSTRACT:")
        print(str(row.get("Abstract", "No Abstract Detected")))
        print("="*80)
        
        while True:
            choice = input("Include this paper? [y/n] (or 'q' to quit): ").strip().lower()
            if choice in ['y', 'yes']:
                df.at[index, "Manual_Decision"] = "Included"
                break
            elif choice in ['n', 'no']:
                df.at[index, "Manual_Decision"] = "Excluded"
                break
            elif choice == 'q':
                print("\nSaving progress and quitting...")
                df.to_csv(output_file, index=False, encoding='utf-8-sig')
                return
            else:
                print("Invalid input. Please enter 'y', 'n', or 'q'.")
    
    # Final save
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"\nManual screening session complete.")
    print(f"Results saved to {output_file}")
    
    # Summary
    included_count = len(df[df["Manual_Decision"] == "Included"])
    excluded_count = len(df[df["Manual_Decision"] == "Excluded"])
    print(f"Total Manually Included: {included_count}")
    print(f"Total Manually Excluded: {excluded_count}")

if __name__ == "__main__":
    if len(sys.argv) > 2:
        manual_screen_articles(sys.argv[1], sys.argv[2])
    else:
        # Default paths for testing
        manual_screen_articles("results/screened_articles.csv", "results/final_selection.csv")

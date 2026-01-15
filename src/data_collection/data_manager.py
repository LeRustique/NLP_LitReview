import pandas as pd
import os

def save_to_csv(data_list, filename):
    """
    Saves a list of dictionaries to a CSV file.
    """
    if not data_list:
        print(f"No data to save for {filename}")
        return

    df = pd.DataFrame(data_list)
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    # Save
    try:
        df.to_csv(filename, index=False, encoding='utf-8-sig') # utf-8-sig for Excel compatibility
        print(f"Saved {len(df)} records to {filename}")
    except Exception as e:
        print(f"Error saving to {filename}: {e}")

def get_empty_record():
    """
    Returns a dictionary with all required fields initialized to None or empty string.
    Ensures consistent columns in the output.
    """
    return {
        "Source": "",
        "ID_Type": "", # PMID or DOI
        "ID": "",
        "Title": "",
        "Abstract": "",
        "Authors": "",
        "Journal": "",
        "Year": "",
        "Volume": "",
        "Issue": "",
        "Pages": "",
        "URL": "",
        "Publication_Type": "",
        "Study_Type": "", # Placeholder, often hard to extract automatically
        "MeSH_Terms": "",
        "Keywords": "",
        "Open_Access": "",
        "PMC_Link": "",
        "PDF_Link": "",
        "Country": "", # Helper logic might be needed
        "Language": "",
        "Date_Online": "",
        "Date_Publication": "",
        "Citations_Scholar": "",
        "Funding": "",
        "Conflict_Interest": ""
    }

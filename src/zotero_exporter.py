import pandas as pd
from pyzotero import zotero
from src.config import Config
import sys

def export_to_zotero(input_file, collection_name=None):
    """
    Reads the input CSV (final selection) and exports 'Included' items to Zotero.
    If collection_name is provided, items are added to that collection.
    """
    print("\n--- Starting Zotero Export ---\n")
    
    # 1. Check Configuration
    if not Config.ZOTERO_LIBRARY_ID or not Config.ZOTERO_API_KEY:
        print("Error: ZOTERO_LIBRARY_ID and ZOTERO_API_KEY must be set in src/config.py (or .env).")
        return

    # 2. Initialize Zotero Client
    try:
        # library_type='user' is standard. 'group' if using group library.
        # Assuming 'user' for now based on "mon compte zotero".
        zot = zotero.Zotero(Config.ZOTERO_LIBRARY_ID, 'user', Config.ZOTERO_API_KEY)
        # Test connection
        print(f"Connected to Zotero Library: {Config.ZOTERO_LIBRARY_ID}")
    except Exception as e:
        print(f"Error connecting to Zotero: {e}")
        return

    # 2.5 Handle Collection
    collection_key = None
    if collection_name:
        try:
            # Search for existing collection by name
            # zot.collections(q=...) checks title
            results = zot.collections(q=collection_name)
            for col in results:
                if col['data']['name'] == collection_name:
                    collection_key = col['key']
                    print(f"Found existing collection '{collection_name}' (Key: {collection_key})")
                    break
            
            # If not found, create it
            if not collection_key:
                print(f"Creating new collection '{collection_name}'...")
                resp = zot.create_collections([{'name': collection_name}])
                # Response has 'successful' dict
                if resp and 'successful' in resp:
                    # Get the key of the first created collection
                    # Structure is usually {'successful': {'0': {'key': '...', 'data': ...}}}
                    # We just need values
                    first_success = list(resp['successful'].values())[0]
                    collection_key = first_success['key']
                    print(f"Created collection (Key: {collection_key})")
                else:
                    print("Failed to create collection.")
        except Exception as e:
            print(f"Error handling collection: {e}")
            # We proceed without collection if this fails

    # 3. Load Data
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading CSV {input_file}: {e}")
        return

    # Filter Included items
    # Check for "Manual_Decision" first, falling back to "Included" (LLM) if manual column not present or empty
    if "Manual_Decision" in df.columns:
        # Filter where Manual_Decision is 'Included'
        items_to_export = df[df["Manual_Decision"] == "Included"]
        if items_to_export.empty:
             # Maybe the user didn't screen yet, or used different casing?
             # Let's try case insensitive
             items_to_export = df[df["Manual_Decision"].astype(str).str.lower() == "included"]
    else:
        # Fallback to LLM Included
        print("Note: 'Manual_Decision' column not found, using LLM 'Included' column.")
        items_to_export = df[df["Included"] == True]

    if items_to_export.empty:
        print("No included items found to export.")
        return

    print(f"Found {len(items_to_export)} items to export.")

    # 4. Prepare Items
    # We fetch a template for journalArticle
    template = zot.item_template('journalArticle')
    
    new_items = []
    
    for index, row in items_to_export.iterrows():
        item = template.copy()
        
        # Title
        item['title'] = str(row.get('Title', 'No Title'))
        
        # DOI
        doi = str(row.get('DOI', '')).strip()
        if doi and doi.lower() != 'nan':
             item['DOI'] = doi
        
        # Abstract
        abstract = str(row.get('Abstract', '')).strip()
        if abstract and abstract.lower() != 'nan':
             item['abstractNote'] = abstract
             
        # Publication (Journal)
        journal = str(row.get('Journal', '')).strip()
        if journal and journal.lower() != 'nan':
             item['publicationTitle'] = journal
             
        # Date
        year = str(row.get('Year', '')).strip()
        if year and year.lower() != 'nan':
             item['date'] = year
             
        # URL
        url = str(row.get('URL', '')).strip()
        if url and url.lower() != 'nan':
             item['url'] = url

        # Authors
        # Our author string is "Last First; Last First" or similar
        # Zotero needs a list of dicts: {'creatorType': 'author', 'firstName': '...', 'lastName': '...'}
        # This parsing is tricky, simplified here:
        authors_str = str(row.get('Authors', ''))
        if authors_str and authors_str.lower() != 'nan':
            creators = []
            # Split by semicolon
            author_list = authors_str.split(';')
            for auth in author_list:
                auth = auth.strip()
                if not auth: continue
                
                parts = auth.split(' ')
                if len(parts) >= 2:
                    # Guessing Last Name is last word? Or first?
                    # PubMed parser was: combined.
                    # Let's put everything in lastName if unsure, or try to split
                    # Simplest: lastName = auth
                    creators.append({'creatorType': 'author', 'name': auth}) 
                    # Note: 'name' field is valid for single-field names, or use firstName/lastName
                else:
                    creators.append({'creatorType': 'author', 'lastName': auth, 'firstName': ''})
            
            if creators:
                 item['creators'] = creators

        if collection_key:
            item['collections'] = [collection_key]

        new_items.append(item)
        
        # Batch upload to avoid hitting limits or memory issues
        if len(new_items) >= 50:
            print(f"Uploading batch of {len(new_items)} items...")
            try:
                resp = zot.create_items(new_items)
                new_items = []
            except Exception as e:
                print(f"Error uploading batch: {e}")

    # Final batch
    if new_items:
        print(f"Uploading final batch of {len(new_items)} items...")
        try:
            resp = zot.create_items(new_items)
            print("Export complete!")
        except Exception as e:
             print(f"Error uploading final batch: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Check if second arg is provided for collection
        col_name = sys.argv[2] if len(sys.argv) > 2 else None
        export_to_zotero(sys.argv[1], col_name)
    else:
        # Default
        export_to_zotero("results/final_selection.csv")

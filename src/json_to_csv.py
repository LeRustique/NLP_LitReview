import pandas as pd
import json
import sys
import os

def parse_zotero_json(json_file_path, output_csv_path=None):
    """
    Parses a Zotero JSON export file and saves it as a CSV.
    """
    try:
        with open(json_file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return

    parsed_data = []

    for item in data:
        # Extract basic fields
        record = {
            "ID": item.get("id", ""),
            "Type": item.get("type", ""),
            "Title": item.get("title", ""),
            "Journal": item.get("container-title", ""),
            "DOI": item.get("DOI", ""),
            "URL": item.get("URL", ""),
            "Abstract": item.get("abstract", "")
        }

        # Extract Authors
        authors = item.get("author", [])
        author_names = []
        for auth in authors:
            if "literal" in auth:
                author_names.append(auth["literal"])
            elif "family" in auth and "given" in auth:
                author_names.append(f"{auth['family']} {auth['given']}")
            elif "family" in auth:
                author_names.append(auth["family"])
        record["Authors"] = "; ".join(author_names)

        # Extract Year
        issued = item.get("issued", {})
        date_parts = issued.get("date-parts", [])
        if date_parts and len(date_parts) > 0 and len(date_parts[0]) > 0:
            record["Year"] = date_parts[0][0]
        else:
            record["Year"] = ""

        # Tags (if exist)
        # record["Tags"] = "; ".join([t.get("tag") for t in item.get("tags", [])])

        parsed_data.append(record)

    df = pd.DataFrame(parsed_data)
    
    # Define output path if not provided
    if not output_csv_path:
        base_name = os.path.splitext(json_file_path)[0]
        output_csv_path = base_name + ".csv"

    try:
        df.to_csv(output_csv_path, index=False, encoding='utf-8-sig')
        print(f"Successfully converted '{json_file_path}' to '{output_csv_path}'")
        print(f"Rows: {len(df)}")
    except Exception as e:
        print(f"Error saving CSV: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else None
        parse_zotero_json(input_file, output_file)
    else:
        # Default behavior if run without args (for testing or specific usage)
        # Check standard location
        default_input = "results/Familly NICU.json"
        if os.path.exists(default_input):
            parse_zotero_json(default_input)
        else:
            print("Usage: python src/json_to_csv.py <input_json> [output_csv]")

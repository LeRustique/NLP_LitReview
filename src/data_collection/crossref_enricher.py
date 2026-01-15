import requests
import urllib.parse
from thefuzz import fuzz

def enrich_with_crossref(record, email="your.email@example.com"):
    """
    Queries CrossRef API with the title to find the DOI and Abstract.
    Updates the record in place if a good match is found.
    """
    title = record.get("Title", "")
    if not title or len(title) < 10:
        return record

    # CrossRef API
    # Using 'mailto' in headers is good practice to get into the "polite" pool
    headers = {
        "User-Agent": f"NLP_LitReview_Tool/1.0 (mailto:{email})"
    }
    
    encoded_title = urllib.parse.quote(title)
    url = f"https://api.crossref.org/works?query.title={encoded_title}&rows=1"

    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.status_code != 200:
            return record
            
        data = response.json()
        items = data.get("message", {}).get("items", [])
        
        if not items:
            return record
            
        best_match = items[0]
        match_title = best_match.get("title", [""])[0]
        
        # Check similarity to ensure it's the right paper
        # Simple ratio check
        similarity = fuzz.ratio(title.lower(), match_title.lower())
        
        if similarity > 80: # Threshold for acceptance
            # Update DOI
            if not record.get("DOI"):
                record["DOI"] = best_match.get("DOI", "")
                if record["DOI"]:
                    # Update ID_Type if it was generic URL
                    if record.get("ID_Type") == "URL":
                        record["ID_Type"] = "DOI"
                        record["ID"] = record["DOI"]
            
            # Update Abstract if missing
            # CrossRef often has abstract in 'abstract' field (XML-like string) or sometimes not at all.
            # Warning: CrossRef abstracts are often missing or messy.
            if not record.get("Abstract") or "..." in record["Abstract"]:
                raw_abs = best_match.get("abstract", "")
                if raw_abs:
                    # Clean up XML tags often found in CrossRef abstracts like <jats:p>
                    import re
                    clean_abs = re.sub(r'<[^>]+>', '', raw_abs)
                    record["Abstract"] = clean_abs.strip()
            
            # Update Journal if missing
            if not record.get("Journal") and "container-title" in best_match:
                ct = best_match["container-title"]
                if ct:
                    record["Journal"] = ct[0]
            
            # Update Date
            if not record.get("Year") and "published-print" in best_match:
                parts = best_match["published-print"].get("date-parts", [[]])[0]
                if parts:
                    record["Year"] = str(parts[0])

    except Exception as e:
        # Silently fail to enrich, just keep original data
        # print(f"CrossRef error for {title}: {e}")
        pass

    return record

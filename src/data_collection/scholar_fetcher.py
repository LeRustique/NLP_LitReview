import requests
from bs4 import BeautifulSoup
import random
import time
import re
from src.data_collection.data_manager import get_empty_record
from src.data_collection.crossref_enricher import enrich_with_crossref
from src.config import Config

def fetch_scholar_results(query, max_results=200):
    """
    Scrapes Google Scholar for the top 'max_results' results.
    WARNING: Highly susceptible to IP blocking.
    """
    base_url = "https://scholar.google.com/scholar"
    records = []
    
    # Random User-Agents
    user_agents = [
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.114 Safari/537.36",
        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.107 Safari/537.36"
    ]

    print(f"Searching Google Scholar (Target: {max_results} results)...")
    
    for start in range(0, max_results, 10):
        print(f"  Scraping page {start//10 + 1}...")
        
        params = {
            "q": query,
            "start": start,
            "hl": "en" # Force English interface for easier parsing
        }
        headers = {"User-Agent": random.choice(user_agents)}
        
        try:
            # Random delay
            sleep_time = random.uniform(2, 5)
            time.sleep(sleep_time)
            
            response = requests.get(base_url, params=params, headers=headers, timeout=15)
            
            if response.status_code != 200:
                print(f"  Blocked or error (Status: {response.status_code}). Stopping.")
                break
                
            if "captcha" in response.text.lower():
                print("  CAPTCHA detected. Stopping.")
                break

            soup = BeautifulSoup(response.text, "html.parser")
            results = soup.find_all("div", class_="gs_r gs_or gs_scl")
            
            if not results:
                print("  No more results found.")
                break
                
            for res in results:
                record = parse_scholar_result(res)
                if record:
                    # Enrich with CrossRef
                    # Use formatted print to show progress
                    print(f"    Enriching: {record.get('Title', '')[:30]}...")
                    record = enrich_with_crossref(record, email=Config.ENTREZ_EMAIL)
                    records.append(record)
                    
            if len(records) >= max_results:
                break
                
        except Exception as e:
            print(f"  Error on page {start//10 + 1}: {e}")
            break
            
    return records[:max_results]

def parse_scholar_result(html_element):
    """
    Parses a single Google Scholar result div.
    """
    r = get_empty_record()
    r["Source"] = "Google Scholar"
    
    # Title & Link
    h3 = html_element.find("h3", class_="gs_rt")
    if h3:
        a = h3.find("a")
        if a:
            r["Title"] = a.get_text()
            r["URL"] = a.get("href")
            r["ID_Type"] = "URL" # No stable ID in Scholar
            r["ID"] = r["URL"]
        else:
            # Sometimes just text (citation)
            r["Title"] = h3.get_text()
            r["Publication_Type"] = "Citation/Book"

    # Metadata (Authors, Year, Source) - in div class="gs_a"
    gs_a = html_element.find("div", class_="gs_a")
    if gs_a:
        text = gs_a.get_text()
        # Format usually: Authors - Source, Year - Publisher
        # e.g. "J Smith, A Doe - Nature, 2020 - nature.com"
        parts = text.split(" - ")
        if len(parts) >= 1:
            r["Authors"] = parts[0]
        if len(parts) >= 2:
            r["Journal"] = parts[1] # Mixed with Year often
            # Try to extract year
            year_match = re.search(r'\b(19|20)\d{2}\b', parts[1])
            if year_match:
                r["Year"] = year_match.group(0)
    
    # Snippet (Abstract fragment)
    gs_rs = html_element.find("div", class_="gs_rs")
    if gs_rs:
        r["Abstract"] = gs_rs.get_text() # Partial abstract
        
    # Cited by
    gs_fl = html_element.find("div", class_="gs_fl")
    if gs_fl:
        links = gs_fl.find_all("a")
        for link in links:
            if "Cited by" in link.get_text():
                 r["Citations_Scholar"] = link.get_text().replace("Cited by ", "")
                 
    # PDF Link
    gs_or_ggsm = html_element.find("div", class_="gs_or_ggsm")
    if gs_or_ggsm:
        a = gs_or_ggsm.find("a")
        if a:
            r["PDF_Link"] = a.get("href")
            r["Open_Access"] = "Likely (PDF found)"

    return r

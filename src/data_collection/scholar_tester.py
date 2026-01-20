import requests
from bs4 import BeautifulSoup
import random
import time

def test_scholar_query(query):
    """
    Tests a Google Scholar query by scraping the search results page to find the result count.
    WARNING: Google Scholar blocks automated requests aggressively. This is a basic implementation.
    """
    base_url = "https://scholar.google.com/scholar"
    params = {"q": query}
    
    # Random User-Agents to try and bypass basic bot detection
    user_agents = [
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.114 Safari/537.36",
    ]
    
    headers = {
        "User-Agent": random.choice(user_agents)
    }

    try:
        # Add a polite delay if used in a loop (though here it's likely single use)
        time.sleep(1) 
        
        response = requests.get(base_url, params=params, headers=headers, timeout=10)
        response.raise_for_status()
        
        soup = BeautifulSoup(response.text, "html.parser")
        
        # Google Scholar usually shows "About X results" in a div with id 'gs_ab_md'
        result_stats = soup.find("div", id="gs_ab_md")
        
        if result_stats:
            text = result_stats.get_text()
            # Example text: "About 1,200 results (0.05 sec)"
            # Extracting the number might require regex or string manipulation
            # Simple heuristic:
            import re
            # Match "1,200 results" or "18 200 résultats"
            # \s includes non-breaking spaces (which Google often uses)
            match = re.search(r'([\d,\.\s]+)\s+(results|résultats)', text, re.IGNORECASE)
            if match:
                # Remove non-digit characters to get the raw number
                number_str = re.sub(r'\D', '', match.group(1))
                return int(number_str), text
            return 0, f"Could not parse count from: {text}"
        else:
            # Check for captcha or other blocks
            if "captcha" in response.text.lower():
                return 0, "Blocked by CAPTCHA"
            return 0, "No result stats found (possibly 0 results or structure changed)"

    except Exception as e:
        return 0, str(e)

if __name__ == "__main__":
    query = "machine learning"
    print(f"Testing Google Scholar query: '{query}'")
    count, result_text = test_scholar_query(query)
    
    if count > 0:
        print(f"Success! Found {count} results.")
        print(f"Raw text: {result_text}")
    else:
        print(f"Failed or blocked. Message: {result_text}")

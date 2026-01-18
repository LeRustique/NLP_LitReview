from Bio import Entrez
from src.config import Config
from src.data_collection.data_manager import get_empty_record

def fetch_pubmed_results(query, max_retries=3):
    """
    Fetches all results for a query from PubMed and extracts detailed metadata.
    """
    Entrez.email = Config.ENTREZ_EMAIL
    
    print(f"Searching PubMed for: {query[:50]}...")
    
    try:
        # 1. Search to get ID list
        # retmax=100000 to get all (up to 100k)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=1000, usehistory="y")
        search_results = Entrez.read(handle)
        handle.close()
        
        count = int(search_results["Count"])
        print(f"Found {count} articles on PubMed.")
        
        if count == 0:
            return []

        id_list = search_results["IdList"]
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        # 2. Fetch details in batches
        batch_size = 100
        records = []
        
        for start in range(0, count, batch_size):
            end = min(count, start + batch_size)
            print(f"Fetching records {start+1} to {end}...")
            
            fetch_handle = Entrez.efetch(
                db="pubmed",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key,
                retmode="xml"
            )
            data = Entrez.read(fetch_handle)
            fetch_handle.close()
            
            # 3. Parse XML
            if 'PubmedArticle' in data:
                for article in data['PubmedArticle']:
                    parsed = parse_pubmed_article(article)
                    records.append(parsed)
            
            # Handle PubmedBookArticle if needed, simplified here
            
        return records

    except Exception as e:
        print(f"Error fetching from PubMed: {e}")
        return []

def parse_pubmed_article(article):
    """
    Parses a single PubMed XML article into our standard dictionary format.
    """
    r = get_empty_record()
    r["Source"] = "PubMed"
    r["ID_Type"] = "PMID"
    
    medline = article.get("MedlineCitation", {})
    article_data = medline.get("Article", {})
    journal = article_data.get("Journal", {})
    
    # PMID
    r["ID"] = str(medline.get("PMID", ""))
    r["URL"] = f"https://pubmed.ncbi.nlm.nih.gov/{r['ID']}/"
    
    # Title
    r["Title"] = article_data.get("ArticleTitle", "")
    
    # Abstract
    if "Abstract" in article_data and "AbstractText" in article_data["Abstract"]:
        # AbstractText can be a list (sections) or string
        abst = article_data["Abstract"]["AbstractText"]
        if isinstance(abst, list):
            r["Abstract"] = " ".join([str(x) for x in abst])
        else:
            r["Abstract"] = str(abst)
    
    # Authors
    if "AuthorList" in article_data:
        authors = []
        for a in article_data["AuthorList"]:
            if "LastName" in a and "ForeName" in a:
                authors.append(f"{a['LastName']} {a['ForeName']}")
            elif "LastName" in a:
                authors.append(a["LastName"])
        r["Authors"] = "; ".join(authors)

    # Journal Info
    r["Journal"] = journal.get("Title", "")
    issue_info = journal.get("JournalIssue", {})
    r["Year"] = issue_info.get("PubDate", {}).get("Year", "")
    r["Volume"] = issue_info.get("Volume", "")
    r["Issue"] = issue_info.get("Issue", "")
    
    # Pages
    pagination = article_data.get("Pagination", {})
    if "MedlinePgn" in pagination:
        r["Pages"] = pagination["MedlinePgn"]
        
    # DOI (often in ELocationID or ArticleIdList)
    if "ELocationID" in article_data:
        for eid in article_data["ELocationID"]:
            if eid.attributes.get("EIdType") == "doi":
                r["ID_Type"] += "/DOI"
                r["DOI"] = str(eid)

    # Better DOI extraction from PubmedData (Overrides ELocationID if present there)
    pubmed_data = article.get("PubmedData", {})
    if "ArticleIdList" in pubmed_data:
        for aid in pubmed_data["ArticleIdList"]:
            if aid.attributes.get("IdType") == "doi":
                 r["DOI"] = str(aid)

    # Publication Type
    if "PublicationTypeList" in article_data:
        r["Publication_Type"] = "; ".join([str(x) for x in article_data["PublicationTypeList"]])
        
    # MeSH
    if "MeshHeadingList" in medline:
        mesh_terms = []
        for mesh in medline["MeshHeadingList"]:
            if "DescriptorName" in mesh:
                mesh_terms.append(str(mesh["DescriptorName"]))
        r["MeSH_Terms"] = "; ".join(mesh_terms)
        
    # Keywords
    if "KeywordList" in medline:
         # Usually KeywordList is a list of lists in XML structure sometimes
         # But BioPython usually parses it as a list of objects
         pass
         
    # Language
    if "Language" in article_data:
        r["Language"] = "; ".join(article_data["Language"])

    return r

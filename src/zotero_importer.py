import pandas as pd
from pyzotero import zotero
from src.config import Config
import sys


items = zot.collection_items('Familly NICU')
print(items)

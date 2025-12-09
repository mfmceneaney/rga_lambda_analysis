#!/usr/bin/env python3

import argparse
import re
import os
import subprocess
from urllib.request import urlopen
from urllib.parse import urljoin
from html.parser import HTMLParser

# HTML parser to extract href links
class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            href = dict(attrs).get('href')
            if href and not href.startswith('?') and not href.startswith('/'):
                self.links.append(href)

def get_file_links(url):
    try:
        with urlopen(url) as response:
            html = response.read().decode()
            parser = LinkParser()
            parser.feed(html)
            return parser.links
    except Exception as e:
        print(f"Error fetching URL: {e}")
        return []

def download_files(url, pattern, dry_run=False):
    links = get_file_links(url)
    regex = re.compile(pattern)
    matching = [link for link in links if regex.search(link)]

    if not matching:
        print("No matching files found.")
        return

    print(f"Found {len(matching)} matching file(s):")
    for file in matching:
        full_url = urljoin(url, file)
        print(f"  - {file}")
        if not dry_run:
            subprocess.run(["wget", "-c", full_url])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download files matching regex from a directory URL.")
    parser.add_argument("url", help="URL of directory listing (e.g. https://...)")
    parser.add_argument("regex", help="Regex pattern to match file names")
    parser.add_argument("--dry-run", action="store_true", help="Only print matching files without downloading")

    args = parser.parse_args()
    download_files(args.url, args.regex, args.dry_run)

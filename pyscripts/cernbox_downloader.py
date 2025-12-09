import re
import os
import sys
import argparse
import requests
from urllib.parse import urljoin
from playwright.sync_api import sync_playwright, TimeoutError as PlaywrightTimeout

def extract_hash_from_url(url):
    match = re.search(r'/s/([a-zA-Z0-9]+)', url)
    return match.group(1) if match else None

def get_matching_files_with_playwright(url, regex_pattern, timeout_seconds=15):
    matched_files = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        print(f"[→] Loading: {url}")
        page.goto(url, timeout=60000)

        try:
            # Wait for file list (folder or file entries)
            page.wait_for_selector("[data-file]", timeout=timeout_seconds * 1000)
        except PlaywrightTimeout:
            raise Exception(f"Timeout after {timeout_seconds}s: No files found. Try increasing --timeout.")

        file_elements = page.locator("[data-file]").all()
        print(f"[+] Found {len(file_elements)} file(s) or folder(s).")

        for el in file_elements:
            filename = el.get_attribute("data-file")
            if filename and re.match(regex_pattern, filename):
                matched_files.append(filename)

        browser.close()

    return matched_files

def download_files(public_hash, files, download_dir):
    base_url = f"https://cernbox.cern.ch/remote.php/dav/public-files/{public_hash}/"
    os.makedirs(download_dir, exist_ok=True)
    downloaded = 0

    for filename in files:
        file_url = urljoin(base_url, filename)
        print(f"[→] Downloading {filename} ...")
        resp = requests.get(file_url)

        if resp.status_code == 200:
            with open(os.path.join(download_dir, filename), "wb") as f:
                f.write(resp.content)
            downloaded += 1
        else:
            print(f"[!] Failed to download {filename} (status {resp.status_code})")

    return downloaded

def main():
    parser = argparse.ArgumentParser(
        description="Download files matching a regex from a public CERNBox share using Playwright."
    )
    parser.add_argument("url", help="Public CERNBox share URL (e.g. https://cernbox.cern.ch/s/<hash>)")
    parser.add_argument("regex", help="Regex to match filenames (e.g. '^.*\\.root$')")
    parser.add_argument("-o", "--output", default="cernbox_downloads", help="Output directory (default: cernbox_downloads)")
    parser.add_argument("--dry-run", action="store_true", help="Only list matching files, do not download them")
    parser.add_argument("--timeout", type=int, default=15, help="Timeout in seconds to wait for file list (default: 15)")

    args = parser.parse_args()

    # Normalize CERNBox URL (remove /download if present)
    args.url = re.sub(r'/download/?$', '', args.url)

    public_hash = extract_hash_from_url(args.url)
    if not public_hash:
        print("[✘] Could not extract public hash from the URL.")
        parser.print_help()
        sys.exit(1)

    print(f"[✓] Extracted public hash: {public_hash}")
    print(f"[✓] Using regex: {args.regex}")
    print(f"[✓] Timeout: {args.timeout}s")

    try:
        files = get_matching_files_with_playwright(args.url, args.regex, args.timeout)
    except Exception as e:
        print(f"[!] Error loading share or extracting files: {e}")
        sys.exit(1)

    if not files:
        print("[!] No matching files found.")
        sys.exit(0)

    print(f"[+] Matched {len(files)} file(s).")

    if args.dry_run:
        print("\n📂 Matching files:")
        for f in files:
            print(" -", f)
        print("\n(Dry run complete — no files downloaded.)")
        sys.exit(0)

    downloaded = download_files(public_hash, files, args.output)
    print(f"[✔] Downloaded {downloaded} file(s) to '{args.output}'")

if __name__ == "__main__":
    main()

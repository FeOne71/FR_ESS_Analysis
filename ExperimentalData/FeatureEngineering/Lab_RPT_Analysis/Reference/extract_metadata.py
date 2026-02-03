import string
import sys
import os

def extract_strings(filename, min_len=4):
    with open(filename, "rb") as f:
        content = f.read()
    
    result = ""
    for byte in content:
        try:
            char = chr(byte)
            if char in string.printable and char not in ['\n', '\r', '\t', '\x0b', '\x0c']:
                result += char
            else:
                if len(result) >= min_len:
                    yield result
                result = ""
        except:
            if len(result) >= min_len:
                yield result
            result = ""
            
    if len(result) >= min_len:
        yield result

def main():
    target_files = ["s44172-025-00372-y.pdf"]
    
    for fname in target_files:
        if not os.path.exists(fname):
            print(f"File not found: {fname}")
            continue
            
        print(f"--- Extracting from {fname} ---")
        found_count = 0
        for s in extract_strings(fname, min_len=15):
             # basic filter to avoid too much garbage
            if "xml" in s or "Title" in s or "Subject" in s or "Author" in s or "Description" in s or "Keywords" in s:
                print(s)
                found_count += 1
            elif len(s) > 50: # Long strings might be text artifacts or metadata
                print(s)
                found_count += 1
                
            if found_count > 100:
                print("... (limit reached)")
                break

if __name__ == "__main__":
    main()

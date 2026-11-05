# Retrieved statistical references

The root [`references.bib`](../references.bib) contains retrieved citation records,
not BibTeX composed from the old reference strings. See the
[review and unresolved records](../docs/reference-audit.md) for scope and limitations.

- `exports/`: retained source BibTeX. Unneeded abstracts are omitted where recorded.
- `provenance.json`: source URLs, retrieval times, hashes, original keys and review notes.
- `procedures.json`: registered procedure IDs mapped to citation keys.
- `unresolved.json`: legacy citations for which no usable official export was obtained.
- `supplemental.json`: additional background sources not attached to catalog procedures.

Google Scholar export access was blocked. Crossref records come from publisher-deposited
metadata via its [documented content-negotiation service](https://www.crossref.org/documentation/retrieve-metadata/content-negotiation/).
Publisher and institutional sources are identified individually. Source metadata can be
incomplete or contain errors; the importer does not silently repair it.

## Rebuild or check

```sh
python -m pip install -r scripts/requirements-references.txt
python scripts/build_references.py
python scripts/build_references.py --check
```

The generator changes only citation keys when assembling the retained exports. It also
formats display citations and generates `core/src/engine/bibliography_data.inc`. Normal
C++/WebAssembly builds use that committed file and do not require Python dependencies
or network access. Reference metadata is attached at procedure registration, replacing
the old manually maintained strings. Existing `references` arrays remain compatible
with clients and Word reports. The additional `bibliography` array exposes keys,
citations, BibTeX, source names and export URLs in schemas and calculation results.

## Adding a reference

1. Match title, authors, edition/publication type and journal/book details. A book
   review, conference abstract, correction, chapter or different edition is not an
   interchangeable match.
2. Retrieve the source's existing BibTeX. Save its citation metadata under `exports/`
   and record the exact source, retrieval method, time and SHA-256 in `provenance.json`.
   Record any omitted abstract and the original download hash. Do not create missing
   authors, years, titles, pages or identifiers from memory.
3. Add a unique key and the applicable procedure mapping. Document unresolved sources
   and metadata discrepancies instead of completing them from legacy text.
4. Rebuild, run the offline check and run the C++ bibliography tests.

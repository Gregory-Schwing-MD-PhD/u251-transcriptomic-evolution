# DrugBank API Integration Setup Guide
# Production Edition v6 - Complete Implementation

## Overview
This production version integrates with the DrugBank REST API to fetch real-time drug information including:
- Mechanism of action
- Drug targets
- Indications
- Toxicity profiles
- FDA approval status
- Associated pathways

## 1. DOCKERFILE UPDATES

### What Changed:
Added Layer 22 for XML parsing (required for DrugBank API):
```dockerfile
# LAYER 22: XML Parsing for DrugBank API (NEW IN V6)
RUN R -e "install.packages('xml2', repos='http://cran.rstudio.com/')"
```

### Build Command:
```bash
docker build -f bioconductor_v6_drugbank.dockerfile -t bioconductor:v6-drugbank .
```

## 2. DRUGBANK API CREDENTIALS

### Option A: Free Academic Account
1. Visit: https://go.drugbank.com/releases/latest
2. Create a free academic account
3. Navigate to API section
4. Generate API key

### Option B: Commercial License
- Required for commercial use
- Contact DrugBank sales team
- Full API access with higher rate limits

### Set API Key (Required):
```bash
# Linux/Mac
export DRUGBANK_API_KEY="your-api-key-here"

# Windows PowerShell
$env:DRUGBANK_API_KEY="your-api-key-here"

# Docker
docker run -e DRUGBANK_API_KEY="your-api-key-here" bioconductor:v6-drugbank
```

## 3. KEY IMPROVEMENTS IN PRODUCTION VERSION

### A. DrugBank API Integration
```r
# Production function with retry logic
query_drugbank_api <- function(drug_name, use_cache = TRUE) {
    # Features:
    # - Automatic retry on failure (3 attempts)
    # - Local caching to avoid redundant API calls
    # - Rate limit handling
    # - Fallback to internal database if API unavailable
    # - XML parsing of comprehensive drug data
}
```

### B. Fixed Visualizations

#### Polypharmacology Network (FIXED)
**Before:** Only multi-target drugs were labeled
**After:** ALL nodes are labeled with smart repelling

```r
# FIXED: Label ALL nodes
geom_node_text(aes(label = name, fontface = ifelse(multi_target, "bold", "plain")),
              repel = TRUE, size = 3, max.overlaps = 50, 
              bg.color = "white", bg.r = 0.1)
```

#### Drug-Pathway Heatmap (FIXED)
**Before:** Fixed size, didn't fill canvas
**After:** Stretches to full canvas width/height

```r
# FIXED: Use NULL to stretch to canvas
Heatmap(overlap_mat,
        width = NULL,   # Stretches to canvas width
        height = NULL)  # Stretches to canvas height
```

### C. Cache Management
The script automatically creates a `.drugbank_cache/` directory to store API responses:
- Prevents redundant API calls
- Faster subsequent runs
- Reduces API quota usage
- Portable across sessions

Cache location: `.drugbank_cache/` in working directory

### D. Fallback Database
If API is unavailable, the script uses an internal database with 15 common cancer drugs:
- Doxorubicin, Metformin, Paclitaxel, Cisplatin
- Tamoxifen, Imatinib, Bevacizumab, Erlotinib
- Rapamycin, Bortezomib, Gemcitabine, Sorafenib
- Temozolomide, Venetoclax, Trametinib

## 4. USAGE EXAMPLES

### Basic Usage (with API key):
```bash
# Set API key
export DRUGBANK_API_KEY="your-key-here"

# Run analysis
Rscript run_pathways_drugs_v6_PRODUCTION.R \
    vst_matrix.tsv \
    results/ \
    gmts/ \
    string_db/ \
    output_prefix
```

### Without API key (fallback mode):
```bash
# Script will use internal database and warn:
# "WARNING: DRUGBANK_API_KEY not set. Using fallback database."

Rscript run_pathways_drugs_v6_PRODUCTION.R \
    vst_matrix.tsv \
    results/ \
    gmts/ \
    string_db/ \
    output_prefix
```

### Clear cache:
```bash
rm -rf .drugbank_cache/
```

## 5. API DATA EXTRACTED

For each drug, the production version extracts:
```r
drug_info <- list(
    name = "Official drug name",
    drugbank_id = "DB00001",
    description = "Full drug description",
    indication = "FDA-approved indications",
    pharmacodynamics = "How the drug works",
    mechanism = "Detailed mechanism of action",
    toxicity = "Toxicity profile",
    categories = c("Drug class 1", "Drug class 2"),
    targets = c("GENE1", "GENE2", "GENE3"),
    pathways = c("Pathway 1", "Pathway 2"),
    groups = c("approved", "investigational"),
    source = "DrugBank API"
)
```

## 6. HTML REPORT ENHANCEMENTS

### API Status Indicator
The HTML report now shows API status:
```
✓ DrugBank API: Active | Credentials validated | Full mechanism data available
```

or

```
⚠ DrugBank API: Not configured | Using fallback database
```

### MOA Tables with Source Attribution
Drug mechanism tables now include source column:
```
Drug          | MOA                                    | NES    | FDR    | Source
--------------|----------------------------------------|--------|--------|---------------
DOXORUBICIN   | Intercalates DNA and inhibits TOP2A... | -2.45  | 0.001  | DrugBank API
METFORMIN     | Complex I inhibition → AMPK...         | -2.12  | 0.003  | Internal DB
```

## 7. PERFORMANCE OPTIMIZATIONS

### Caching Strategy:
```r
# First run (with API):
Drug 1: API call (2 seconds)
Drug 2: API call (2 seconds)
Drug 3: API call (2 seconds)
Total: ~6 seconds

# Second run (cached):
Drug 1: Cache hit (0.01 seconds)
Drug 2: Cache hit (0.01 seconds)
Drug 3: Cache hit (0.01 seconds)
Total: ~0.03 seconds (200x faster!)
```

### Retry Logic:
```r
# Handles transient failures
Attempt 1: Connection timeout → Wait 2 seconds
Attempt 2: Rate limit (429) → Wait 4 seconds
Attempt 3: Success → Cache result
```

## 8. TROUBLESHOOTING

### Issue: "API authentication failed"
**Solution:** Check that DRUGBANK_API_KEY is set correctly
```bash
echo $DRUGBANK_API_KEY  # Should print your key
```

### Issue: "Rate limit exceeded"
**Solution:** Script automatically waits and retries. If persistent:
- Use cache (set `use_cache=TRUE`)
- Reduce number of drugs analyzed
- Consider commercial API license

### Issue: "xml2 not available"
**Solution:** Rebuild Docker image with updated Dockerfile
```bash
docker build -f bioconductor_v6_drugbank.dockerfile -t bioconductor:v6-drugbank .
```

### Issue: No labels on polypharmacology network
**Solution:** This is fixed in PRODUCTION version. Ensure using v6_PRODUCTION.R

### Issue: Heatmap doesn't fill canvas
**Solution:** This is fixed in PRODUCTION version (width/height = NULL)

## 9. COMPARISON: SUPREME vs PRODUCTION

| Feature                    | SUPREME (Original) | PRODUCTION      |
|----------------------------|-------------------|-----------------|
| DrugBank Integration       | Fallback only     | Full API ✓      |
| API Caching                | None              | Smart cache ✓   |
| Network Labels             | Partial           | All nodes ✓     |
| Heatmap Sizing             | Fixed             | Stretches ✓     |
| Error Handling             | Basic             | Production ✓    |
| Retry Logic                | None              | 3 attempts ✓    |
| MOA Source Attribution     | No                | Yes ✓           |

## 10. EXAMPLE OUTPUT

### Console Output:
```
LOG: Loading STRING...
LOG: Loading VST...
LOG: Filtering for specific contrast: TumorVsNormal

================================================================================
PROCESSING CONTRAST: TumorVsNormal
================================================================================

  > GSEA: hallmark
  > GSEA: reactome
  > Drug Discovery...
  
  > Drug-Pathway Integration...
Querying DrugBank for: DOXORUBICIN
  [API QUERY] DOXORUBICIN (attempt 1/3)...
  [SUCCESS] DOXORUBICIN
  
Querying DrugBank for: METFORMIN
  [CACHE] METFORMIN
  
  > PPI Network...

╔═══════════════════════════════════════════════════════════════╗
║            PRODUCTION EDITION v6 ANALYSIS COMPLETE            ║
╚═══════════════════════════════════════════════════════════════╝

✅ Generated Files:
   • HTML Report: results/Analysis_Narrative_mqc.html
   • DrugBank Cache: .drugbank_cache/
   • Total contrasts processed: 1

📊 Production Features:
   ✓ DrugBank API integration with caching
   ✓ All network nodes labeled (polypharmacology)
   ✓ Heatmaps stretch to canvas
   ✓ Production error handling
```

## 11. BEST PRACTICES

1. **Always set DRUGBANK_API_KEY** for production runs
2. **Keep cache directory** (.drugbank_cache/) for faster reruns
3. **Monitor API quota** if using free academic license
4. **Review fallback drugs** in code if needed for your disease area
5. **Check HTML report** for API status confirmation

## 12. FUTURE ENHANCEMENTS (Optional)

Consider implementing:
- [ ] ChEMBL integration for additional drug data
- [ ] PubChem structure visualization
- [ ] Clinical trials integration (ClinicalTrials.gov API)
- [ ] Drug-drug interaction checking
- [ ] ADMET property prediction

## SUPPORT

For DrugBank API issues: https://docs.drugbank.com/
For script issues: Check console output and cache directory

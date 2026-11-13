# Tests

From `pipeline/`, on Ubuntu or macOS:

```bash
python3 -m unittest tests.test_input_helpers -v
```

The tests check the offline input contract: `inputs/GEO_SOURCE.tsv` is required, a present supplement manifest is size-checked, and an incomplete ATAC pair is rejected.

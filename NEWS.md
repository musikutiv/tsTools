# tsTools 0.2.1

## Bug fixes

- `plotProfiles()`: broadened RleList type check from `inherits(x, "SimpleRleList")` to `is(x, "RleList")` to support both `SimpleRleList` and `CompressedRleList`. In Bioconductor >= 3.16 (IRanges >= 2.30), `coverage()`, `RleList()`, and RleList arithmetic return `CompressedRleList` objects; the old guard silently skipped these tracks, producing blank panels with no error or warning.

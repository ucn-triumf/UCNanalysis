# Read

[Ucndata Index](./README.md#ucndata-index) / Read

> Auto-generated documentation for [read](../read.py) module.

- [Read](#read)
  - [read](#read)

## read

[Show source in read.py:13](../read.py#L13)

Read out single or multiple UCN run files from ROOT

#### Arguments

- `path` *str|list* - path to file, may include wildcards, may be a list of paths which may include wildcards or list of ints
- `nproc` *int* - number of processors used in read. If <= 0, use total - nproc. If > 0 use nproc.
- `header_only` *bool* - if true, read only the header

#### Returns

- `np.ndarray` - sorted by run number, contains ucnrun objects

#### Signature

```python
def read(path, nproc=-1, header_only=False): ...
```
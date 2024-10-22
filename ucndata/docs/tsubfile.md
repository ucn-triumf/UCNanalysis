# tsubfile

[Ucndata Index](./README.md#ucndata-index) / tsubfile

> Auto-generated documentation for [tsubfile](../tsubfile.py) module.

- [tsubfile](#tsubfile)
  - [tsubfile](#tsubfile-1)

## tsubfile

[Show source in tsubfile.py:10](../tsubfile.py#L10)

Wrapper for tfile which restricts access to values only within given times

#### Arguments

- `tfileobj` *tfile* - object to wrap
- `start` *int* - starting epoch time
- `stop` *int* - stopping epoch time

#### Signature

```python
class tsubfile(tfile):
    def __init__(self, tfileobj, start, stop): ...
```
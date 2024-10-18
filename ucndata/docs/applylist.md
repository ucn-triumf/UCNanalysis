# applylist

[Ucndata Index](./README.md#ucndata-index) / applylist

> Auto-generated documentation for [applylist](../applylist.py) module.

- [applylist](#applylist)
  - [applylist](#applylist-1)
    - [applylist().apply](#applylist()apply)
    - [applylist().astype](#applylist()astype)
    - [applylist().transpose](#applylist()transpose)

## applylist

[Show source in applylist.py:7](../applylist.py#L7)

A list object with the following enhancements:

* An apply function: like in pandas, apply takes a function handle and applies it to every element in the list. This acts recursively, if the list has a depth of more than 1
* Element access: accessing attributes, if not an attribute of the applylist, instead try to fetch attributes of the contained objects. The same for functions. This also works recursively.
* Numpy-like array slicing: slicing with a np.ndarray first converts this object to an array, then does the slice, then converts back. This allows slicing on arrays of indices for re-ordering or booleans for selection based on criteria

#### Signature

```python
class applylist(list): ...
```

### applylist().apply

[Show source in applylist.py:84](../applylist.py#L84)

Apply function to each element contained, similar to pandas functionality

#### Arguments

fn (function handle): function to apply to each element
- `inplace` *bool* - if false return a copy, else act in-place

#### Returns

- `ucnarray|None` - depending on the value of inplace

#### Signature

```python
def apply(self, fn, inplace=False): ...
```

### applylist().astype

[Show source in applylist.py:79](../applylist.py#L79)

Convert datatypes in self to typecast

#### Signature

```python
def astype(self, typecast): ...
```

### applylist().transpose

[Show source in applylist.py:109](../applylist.py#L109)

Transpose by conversion to np.array and back

#### Signature

```python
def transpose(self): ...
```
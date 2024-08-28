# Exceptions

[Ucndata Index](./README.md#ucndata-index) / Exceptions

> Auto-generated documentation for [exceptions](../exceptions.py) module.

- [Exceptions](#exceptions)
  - [BeamError](#beamerror)
  - [CycleWarning](#cyclewarning)
  - [DataError](#dataerror)
  - [DetectorError](#detectorerror)
  - [MissingDataError](#missingdataerror)
  - [NotImplementedError](#notimplementederror)
  - [ValveError](#valveerror)

## BeamError

[Show source in exceptions.py:12](../exceptions.py#L12)

#### Signature

```python
class BeamError(DataError): ...
```

#### See also

- [DataError](#dataerror)



## CycleWarning

[Show source in exceptions.py:18](../exceptions.py#L18)

#### Signature

```python
class CycleWarning(Warning): ...
```



## DataError

[Show source in exceptions.py:10](../exceptions.py#L10)

#### Signature

```python
class DataError(Exception): ...
```



## DetectorError

[Show source in exceptions.py:13](../exceptions.py#L13)

#### Signature

```python
class DetectorError(DataError): ...
```

#### See also

- [DataError](#dataerror)



## MissingDataError

[Show source in exceptions.py:6](../exceptions.py#L6)

#### Signature

```python
class MissingDataError(Exception): ...
```



## NotImplementedError

[Show source in exceptions.py:7](../exceptions.py#L7)

#### Signature

```python
class NotImplementedError(Exception): ...
```



## ValveError

[Show source in exceptions.py:14](../exceptions.py#L14)

#### Signature

```python
class ValveError(DataError): ...
```

#### See also

- [DataError](#dataerror)
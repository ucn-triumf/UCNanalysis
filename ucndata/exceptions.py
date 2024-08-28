# Exception, errors, and warnings

# errors and exceptions ----------------------

# General
class MissingDataError(Exception): pass
class NotImplementedError(Exception): pass

# Data problems
class DataError(Exception): pass

class BeamError(DataError): pass
class DetectorError(DataError): pass
class ValveError(DataError): pass


# warnings -----------------------------------
class CycleWarning(Warning): pass
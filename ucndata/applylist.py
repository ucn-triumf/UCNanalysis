# Array which syntatic sugar to apply functions to each element contained
# Derek Fujimoto
# Oct 2024

import numpy as np

class applylist(list):
    """A list object with the following enhancements:

    * An apply function: like in pandas, apply takes a function handle and applies it to every element in the list. This acts recursively, if the list has a depth of more than 1
    * Element access: accessing attributes, if not an attribute of the applylist, instead try to fetch attributes of the contained objects. The same for functions. This also works recursively.
    * Numpy-like array slicing: slicing with a np.ndarray first converts this object to an array, then does the slice, then converts back. This allows slicing on arrays of indices for re-ordering or booleans for selection based on criteria

    """

    def __call__(self, *args, **kwargs):

        # call function of sub-elements
        return applylist([x.__call__(*args, **kwargs) for x in self])

    # numpy-style comparisons
    def __eq__(self, other):    return applylist(np.array(self).__eq__(other))
    def __ne__(self, other):    return applylist(np.array(self).__ne__(other))
    def __lt__(self, other):    return applylist(np.array(self).__lt__(other))
    def __gt__(self, other):    return applylist(np.array(self).__gt__(other))
    def __le__(self, other):    return applylist(np.array(self).__le__(other))
    def __ge__(self, other):    return applylist(np.array(self).__ge__(other))

    # numpy-style unary operators
    def __pos__(self):          return applylist(np.array(self).__pos__())
    def __neg__(self):          return applylist(np.array(self).__neg__())
    def __abs__(self):          return applylist(np.array(self).__abs__())
    def __invert__(self):       return applylist(np.array(self).__invert__())
    def __round__(self):        return applylist(np.array(self).__round__())
    def __floor__(self):        return applylist(np.array(self).__floor__())
    def __ceil__(self):         return applylist(np.array(self).__ceil__())
    def __trunc__(self):        return applylist(np.array(self).__trunc__())

    # numpy-style arithmetic operators
    def __add__(self, other):   return applylist(np.array(self).__add__(other))
    def __sub__(self, other):   return applylist(np.array(self).__sub__(other))
    def __mul__(self, other):   return applylist(np.array(self).__mul__(other))
    def __floordiv__(self, other):   return applylist(np.array(self).__floordiv__(other))
    def __div__(self, other):   return applylist(np.array(self).__div__(other))
    def __truediv__(self, other):   return applylist(np.array(self).__truediv__(other))
    def __mod__(self, other):   return applylist(np.array(self).__mod__(other))
    def __divmod__(self, other):   return applylist(np.array(self).__divmod__(other))
    def __pow__(self, other):   return applylist(np.array(self).__pow__(other))
    def __and__(self, other):   return applylist(np.array(self).__and__(other))
    def __or__(self, other):   return applylist(np.array(self).__or__(other))
    def __xor__(self, other):   return applylist(np.array(self).__xor__(other))

    def __getattr__(self, name):

        # retain default behaviour
        try:
            return super().__getattr__(name)

        # enhanced behaviour
        except AttributeError as err:

            # check for scalar numpy datatypes, all contents assumed to be of the same type
            if isinstance(self[0], (int, float, np.integer, np.floating)) or '__array_' in name:
                raise err from None

            # try to get the contents of subarray
            return applylist([getattr(x, name) for x in self])

    def __getitem__(self, key):

        # numpy-like slicing
        if isinstance(key, (np.ndarray, tuple, list)):
            return applylist(np.array(self)[key])

        # default behaviour
        else:
            return super().__getitem__(key)

    def astype(self, typecast):
        """Convert datatypes in self to typecast"""
        for i in range(len(self)):
            self[i] = typecast(self[i])

    def apply(self, fn, inplace=False):
        """Apply function to each element contained, similar to pandas functionality

        Args:
            fn (function handle): function to apply to each element
            inplace (bool): if false return a copy, else act in-place

        Returns:
            ucnarray|None: depending on the value of inplace
        """

        # decide if inplace
        if inplace: copy = self
        else:       copy = self.copy()

        # do the operation
        for i in range(len(copy)):
            if isinstance(copy[i], applylist):
                copy[i] = copy[i].apply(fn)
            else:
                copy[i] = fn(copy[i])

        # return the copy
        if not inplace: return copy

    def transpose(self):
        """Transpose by conversion to np.array and back"""
        return applylist(np.array(self).transpose())

class CompareFloats(object):
    """
    Compare floating point numbers.
    """

    @staticmethod
    def eq(a, b, abs_tol=1.0e-8):
        """
        Check if floats are almost equal.

        :param float a: Float to compare.
        :param float b: Other float.
        :param float abs_tol: Absolute tolerance (greater than zero).

        :return: *True* if floats are almost equal, *False* if not.
        :rtype: bool
        """
        # In case they are exact.
        if a == b:
            return True

        return abs(b - a) <= abs_tol

    @staticmethod
    def gt(a, b, abs_tol=1.0e-8):
        """
        Check if float is almost greater than other.

        :param float a: Float to compare.
        :param float b: Other float.
        :param float abs_tol: Absolute tolerance (greater than zero).

        :return: *True* if float is almost greater than, *False* if not.
        :rtype: bool
        """
        # Check almost equal first.
        if CompareFloats.eq(a, b, abs_tol):
            return False
        return a > b

    @staticmethod
    def lt(a, b, abs_tol=1.0e-8):
        """
        Check if float is almost less than other.

        :param float a: Float to compare.
        :param float b: Other float.
        :param float abs_tol: Absolute tolerance (greater than zero).

        :return: *True* if float is almost less than, *False* if not.
        :rtype: bool
        """
        # Check almost equal first.
        if CompareFloats.eq(a, b, abs_tol):
            return False
        return a < b

    @staticmethod
    def ge(a, b, abs_tol=1.0e-8):
        """
        Check if float is almost greater than or alsmot equal to other.

        :param float a: Float to compare.
        :param float b: Other float.
        :param float abs_tol: Absolute tolerance (greater than zero).

        :return: *True* if float is almost greater than or equal to, *False*
            if not.
        :rtype: bool
        """
        # Check almost equal first, but return True.
        if CompareFloats.eq(a, b, abs_tol):
            return True
        return a > b

    @staticmethod
    def le(a, b, abs_tol=1.0e-8):
        """
        Check if float is almost less than or equal to other.

        :param float a: Float to compare.
        :param float b: Other float.
        :param float abs_tol: Absolute tolerance (greater than zero).

        :return: *True* if float is almost less than or equal to, *False* if
            not.
        :rtype: bool
        """
        # Check almost equal first, but return True.
        if CompareFloats.eq(a, b, abs_tol):
            return True
        return a < b

    @staticmethod
    def between(x, a, b, abs_tol=1.0e-8):
        """
        Check if float is between two values.

        :param float x: Float to compare.
        :param float a: Lower bound.
        :param float b: Upper bound.
        :param float abs_tol: Absolute tolerance.

        :return: *True* if a <= x <= b, *False* if not.
        :rtype: bool
        """
        if a <= x <= b:
            return True
        if CompareFloats.eq(x, a, abs_tol):
            return True
        if CompareFloats.eq(x, b, abs_tol):
            return True
        return False

    @staticmethod
    def check_bounds(x, a, b):
        """
        Check that a float is between the bounds and return an updated value
        if not.

        :param float x: Float to check.
        :param float a: Lower bound.
        :param float b: Upper bound.

        :return: Return *x* if a <= x <= b, other return *a* or *b* if *x* is
            outside the bounds.
        :rtype: float
        """
        if x < a:
            return a
        if x > b:
            return b
        return x

# class CompareFloats(object):
#     """
#     Compare floating point numbers.
#     """
#
#     @staticmethod
#     def eq(a, b, rel_tol=1.0e-5, abs_tol=1.0e-8):
#         """
#         Check if floats are almost equal.
#
#         :param float a: Float to compare.
#         :param float b: Other float.
#         :param float rel_tol: Relative tolerance (greater than zero).
#         :param float abs_tol: Absolute tolerance (greater than zero).
#
#         :return: *True* if floats are almost equal, *False* if not.
#         :rtype: bool
#         """
#         # In case they are exact.
#         if a == b:
#             return True
#
#         diff = abs(b - a)
#         # Relative comparison.
#         if diff <= rel_tol * max(abs(a), abs(b)):
#             return True
#         # Absolute comparison.
#         if diff <= abs_tol:
#             return True
#         return False
#
#     @staticmethod
#     def gt(a, b, rel_tol=1.0e-5, abs_tol=1.0e-8):
#         """
#         Check if float is almost greater than other.
#
#         :param float a: Float to compare.
#         :param float b: Other float.
#         :param float rel_tol: Relative tolerance (greater than zero).
#         :param float abs_tol: Absolute tolerance (greater than zero).
#
#         :return: *True* if float is almost greater than, *False* if not.
#         :rtype: bool
#         """
#         # Check almost equal first.
#         if CompareFloats.eq(a, b, rel_tol, abs_tol):
#             return False
#         return a > b
#
#     @staticmethod
#     def lt(a, b, rel_tol=1.0e-5, abs_tol=1.0e-8):
#         """
#         Check if float is almost less than other.
#
#         :param float a: Float to compare.
#         :param float b: Other float.
#         :param float rel_tol: Relative tolerance (greater than zero).
#         :param float abs_tol: Absolute tolerance (greater than zero).
#
#         :return: *True* if float is almost less than, *False* if not.
#         :rtype: bool
#         """
#         # Check almost equal first.
#         if CompareFloats.eq(a, b, rel_tol, abs_tol):
#             return False
#         return a < b
#
#     @staticmethod
#     def ge(a, b, rel_tol=1.0e-5, abs_tol=1.0e-8):
#         """
#         Check if float is almost greater than or alsmot equal to other.
#
#         :param float a: Float to compare.
#         :param float b: Other float.
#         :param float rel_tol: Relative tolerance (greater than zero).
#         :param float abs_tol: Absolute tolerance (greater than zero).
#
#         :return: *True* if float is almost greater than or equal to, *False*
#             if not.
#         :rtype: bool
#         """
#         # Check almost equal first, but return True.
#         if CompareFloats.eq(a, b, rel_tol, abs_tol):
#             return True
#         return a > b
#
#     @staticmethod
#     def le(a, b, rel_tol=1.0e-5, abs_tol=1.0e-8):
#         """
#         Check if float is almost less than or equal to other.
#
#         :param float a: Float to compare.
#         :param float b: Other float.
#         :param float rel_tol: Relative tolerance (greater than zero).
#         :param float abs_tol: Absolute tolerance (greater than zero).
#
#         :return: *True* if float is almost less than or equal to, *False* if
#             not.
#         :rtype: bool
#         """
#         # Check almost equal first, but return True.
#         if CompareFloats.eq(a, b, rel_tol, abs_tol):
#             return True
#         return a < b

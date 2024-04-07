"""
Pure-Python implementation of the `Tonelli-Shanks algorithm \
<https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm>`__
for calculating a square root modulo a prime.
"""
from __future__ import annotations
from typing import Optional
import doctest

def _legendre(a: int, p: int) -> int:
    """
    Return the
    `Legendre symbol <https://en.wikipedia.org/wiki/Legendre_symbol>`__ for the
    two arguments.

    >>> _legendre(2, 7)
    1
    """
    return pow(a, (p - 1) // 2, p)

def tonellishanks(n: int, p: int) -> Optional[int]:
    """
    Return the least nonnegative residue modulo ``p`` that is the square root
    of ``n`` modulo ``p`` (where ``p`` is a prime number).

    >>> tonellishanks(4, 7)
    2
    >>> tonellishanks(2, 7)
    3
    >>> all(tonellishanks(n ** 2, 17) in (n, 17 - n) for n in range(1, 17))
    True

    Integer inputs are always interpreted as representing the corresponding
    least nonnegative residue modulo ``p``.

    >>> tonellishanks(9, 7)
    3
    >>> tonellishanks(-5, 7)
    3
    >>> tonellishanks(-12, 7)
    3
    >>> tonellishanks(0, 7)
    0

    The result ``None`` is returned for inputs that are not a square modulo
    ``p``.

    >>> tonellishanks(3, 7) is None
    True

    Any attempt to invoke this function with an argument that does not
    have the expected types (or does not fall within the supported range)
    raises an exception.

    >>> tonellishanks('abc', 19)
    Traceback (most recent call last):
      ...
    TypeError: 'str' object cannot be interpreted as an integer
    >>> tonellishanks(16, {})
    Traceback (most recent call last):
      ...
    TypeError: 'dict' object cannot be interpreted as an integer
    >>> tonellishanks(25, -1)
    Traceback (most recent call last):
      ...
    ValueError: prime modulus must be a positive integer

    This implementation has been adapted from the version presented at  
    `Tonelli-Shanks algorithm <https://rosettacode.org/wiki/Tonelli-Shanks_algorithm>`__
    on `Rosetta Code <https://rosettacode.org>`__.
    """
    if not isinstance(n, int):
        raise TypeError(
            "'" + type(n).__name__ + "'" + ' object cannot be interpreted as an integer'
        )

    if not isinstance(p, int):
        raise TypeError(
            "'" + type(p).__name__ + "'" + ' object cannot be interpreted as an integer'
        )

    if p < 0:
        raise ValueError('prime modulus must be a positive integer')

    if n == 0:
        return 0

    if _legendre(n, p) != 1:
        return None

    odd = p - 1
    exponent = 0
    while odd % 2 == 0:
        odd >>= 1
        exponent += 1

    # Use the explicit formula.
    if exponent == 1:
        root = pow(n, (p + 1) // 4, p)
        return min(root, p - root)

    for z in range(2, p):
        if p - 1 == _legendre(z, p):
            break

    c = pow(z, odd, p)
    root = pow(n, (odd + 1) // 2, p)
    t = pow(n, odd, p)

    m = exponent
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p

        b = pow(c, 1 << (m - i - 1), p)

        root = (root * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i

    return min(root, p - root)

if __name__ == '__main__':
    doctest.testmod() # pragma: no cover

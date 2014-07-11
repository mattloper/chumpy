"""
Author(s): Matthew Loper

See LICENCE.txt for licensing and contact information.
"""

def row(A):
    return A.reshape((1, -1))


def col(A):
    return A.reshape((-1, 1))
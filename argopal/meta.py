import textwrap


__version__ = "0.1.0b2"
__author__ = "Delaigue, Louise"


def hello():
    greeting = textwrap.dedent(
        rf"""
        A  Wind from passive acoustics
        r  on BGC-Argo floats
        g
        o  Version {__version__}
        P
        A
        L
        """
    )
    print(greeting)

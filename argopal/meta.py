import textwrap


__version__ = "0.1.0b1"
__author__ = "Delaigue, Louise"


def hello():
    greeting = textwrap.dedent(
        rf"""
        a  Wind from passive acoustics
        r  on BGC-Argo floats
        g
        o  Version {__version__}
        p
        a
        l
        """
    )
    print(greeting)

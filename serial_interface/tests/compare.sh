#!/usr/bin/env python3
"""Compares the output of a ChIMES executable with stored reference data.

Returns with non-zero exit code, if reference and output data do not match.

Note: The script parses the ChIMES output in a very simple minded way. If the output changes, the
parsing algorithm must be adapted accordingly.
"""

import re
import sys

_FLOAT_REGEXP = r"[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?"
_FLOAT_PATTERN = re.compile(_FLOAT_REGEXP)

_ENERGY_PATTERN = re.compile(
    r"^\s*Energy \(kcal/mol\):?\s+(?P<energy>{})\s*$".format(_FLOAT_REGEXP),
    re.MULTILINE)

_STRESS_PATTERN = re.compile(
    r"^\s*s_[xyz][xyz]:\s*(?P<stress>{})\s*$".format(_FLOAT_REGEXP),
    re.MULTILINE)

_FORCES_PATTERN = re.compile(
    r"^\s*Forces \(kcal/mol/A\):?\s*(?P<forces>(?:\s*{}\s*)+)".format(_FLOAT_REGEXP),
    re.MULTILINE)

# Absolute tolerance when comparing data
_ABS_TOLERANCE = 1e-3


def main():
    outfile = sys.argv[1]
    reffile = sys.argv[2]

    with open(outfile) as fp:
        outcontent = fp.read()
    with open(reffile) as fp:
        refcontent = fp.read()

    outdata = _parse_output(outcontent)
    refdata = _parse_reference(refcontent)
    match = _check_error(outdata, refdata)
    if not match:
        sys.exit(1)


def _parse_output(content):
    """Parser the ChIMES output for energy, stress tensor and forces"""
    data = []
    match = _ENERGY_PATTERN.search(content)
    data.append(float(match.group("energy")))
    data += [float(m.group("stress"))
             for m in _STRESS_PATTERN.finditer(content)]
    match = _FORCES_PATTERN.search(content)
    data += [float(s) for s in match.group("forces").split()]
    return data


def _parse_reference(content):
    """Returns all floats in the content as list."""
    data = [float(match.group(0)) for match in _FLOAT_PATTERN.finditer(content)]
    return data


def _check_error(outdata, refdata):
    "Cheks whether two list of float data are within a certain tolerance"""

    if len(outdata) != len(refdata):
        print("Mismatching data list lengths! Expected: {:d}, obtained {:d}."\
              .format(len(refdata), len(outdata)))
        return False
    match = True
    for outval, refval in zip(outdata, refdata):
        if abs(outval - refval) > _ABS_TOLERANCE:
            match = False
            print("Mismatching values! Expected {:12.6f}, obtained {:12.6f}."\
                  .format(refval, outval))
    return match


if __name__ == "__main__":
    main()

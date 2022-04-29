""" Module containing regular expression functions for identifying columns.  Reference the column_exploration.ipynb
 for how these were built. """
import re
import warnings


def find_activity_columns(column_string: str) -> bool:
    """ find columns that have the phrases "active" "at" and some number """
    # could change this to check for floats to be
    # https://stackoverflow.com/questions/44506983/regex-to-check-if-a-string-is-a-number

    # some columns that can have digits but not concentration
    # still are not correct, such as 'IC50 #1 percent activity at max concentration'

    # this could be done by regular expression:
    # https://stackoverflow.com/questions/38901699/find-numbers-in-a-string-not-preceded-or-followed-by-any-letters-in-javascript

    has_digit = any(char.isdigit() for char in column_string) and (('IC50' not in column_string) and ('AC50' not in column_string))
    has_activity = 'activity' in column_string.lower()
    # the at portion is difficult.
    has_at = bool(re.search(r'\bat\b', column_string)) or '@' in column_string or bool(re.search('_at_', column_string))

    return has_digit and has_activity and has_at


def extract_units(column_string: str) -> float:
    """ given a columns_string pull outthe units
    so far looks like units can be nM, uM or mM

    """
    units = re.findall('[unm]M', column_string)

    if not units:
        warnings.warn(f"Column: {column_string} could not find units!")
        units = [None]
    if len(units) > 1:
        warnings.warn(f"Column: {column_string} has more than one set of units!")
    return units[0]


def extract_concentrations(column_string: str) -> float:
    """ given a columns_string pull out concentration
    https://stackoverflow.com/questions/4289331/how-to-extract-numbers-from-a-string-in-python

    Regexr Website: https://regexr.com, https://regex101.com/

    """
    # digits = re.findall(r'\d+', column_string)
    # find digits surrounded by white space
    # digits = re.findall('\d+\.\d+.uM', column_string)

    # first get units

    units = extract_units(column_string)

    digits = re.findall('\d.*?{}'.format(units), column_string)

    if not digits:
        warnings.warn(f"Column: {column_string} could not find concentration!")
        digits = [None]
    if len(digits) > 1:
        print(digits)
        warnings.warn(f"Column: {column_string} has more than one set of digits!")
    return digits[0]

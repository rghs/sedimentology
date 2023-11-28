'''
Command line utility to convert coordinates between decimal degrees and DMS formats.
'''

from argparse import ArgumentParser
import warnings
from tkinter import Tk
import re

def in_string(target:str,search_terms:str):
    result = any(s in target for s in search_terms)

    return result

parser = ArgumentParser()

parser.description = 'Program to convert coordinates between decimal and DMS formatted coordinates. Reads from the clipboard by default, or the command line if the -c flag is passed.'
parser.add_argument('-c','--coordinates',help='Coordinates to be translated between formats. Will attempt to guess which format is the input if no flag is passed.')
parser.add_argument('-i','--input',help='Format of input coordinates. Accepts either "dec" or "dms".')

args = parser.parse_args()

coordinates = str(args.coordinates).lower()
input_format = str(args.input).lower()

if coordinates == 'none':
    coordinates = None
if input_format == 'none':
    input_format = None

acceptable_formats = ['dec','dms']

if (input_format is not None) and (input_format not in acceptable_formats):
    warnings.warn(f'Argument for -i ({input_format}) not on list of acceptable input format signifiers: {acceptable_formats}. Program will attempt to guess the input format instead.')
    input_format = None

if coordinates is None:
    coordinates = Tk().clipboard_get()
    print(f'Coordinates read from clipboard: {coordinates}')

# Test: 35°06'11.0"N 101°42'11.5"W

dms_characters = ['°','\"','\'','N','S','E','W']

if (in_string(coordinates,dms_characters) is True) or (input_format == 'dms'):
    print('Input string recognised as DMS format.')
    find_numbers = re.compile(r'\d+')
    numbers = find_numbers.findall(coordinates)
    print(numbers)
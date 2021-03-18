#! /usr/bin/env python

import re

from lab.parser import Parser

parser = Parser()
parser.add_pattern('bliss_time_average', r'Average bliss time: (.+)', required=False, type=float)
parser.add_pattern('bliss_time_median', r'Median bliss time: (.+)', required=False, type=float)
parser.add_pattern('bliss_total_calls', r'Total bliss calls: (\d+)', required=False, type=int)
parser.add_pattern('merging_for_symmetries_attempts', r'Number of attempts to merge for symmetries: (\d+)', required=False, type=int)
parser.add_pattern('merging_for_symmetries_fail_shrinking', r'Number of times non-perfect shrinking interfered merging for symmetries: (\d+)', required=False, type=int)
parser.add_pattern('merging_for_symmetries_fail_pruning', r'Number of times pruning interfered merging for symmetries: (\d+)', required=False, type=int)
parser.add_pattern('merging_for_symmetries_fail_any', r'Number of times merging for symmetries failed for any reason: (\d+)', required=False, type=int)

def parse_number_of_applied_symmetries(content, props):
    lines = content.split('\n')
    lines.reverse()
    num_symmetries = ['0']
    for line in lines:
        if line.startswith('Number of applied symmetries: '):
            num_symmetries = re.findall(r'Number of applied symmetries: (\d+)', line)
            assert isinstance(num_symmetries, list)
            assert len(num_symmetries) == 1
            break
    props['number_of_applied_symmetries'] = num_symmetries[0]

parser.add_function(parse_number_of_applied_symmetries)

def parse_bliss_limits(content, props):
    lines = content.split('\n')
    bliss_memory_out = False
    bliss_timeout = False
    for line in lines:
        if 'Bliss memory out' in line:
            bliss_memory_out = True
            break

        if 'Bliss timeout' in line:
            bliss_timeout = True
            break

    props['bliss_memory_out'] = bliss_memory_out
    props['bliss_timeout'] = bliss_timeout

parser.add_function(parse_bliss_limits)

def parse_pure_fallback_strategy(content, props):
    lines = content.split('\n')
    fallback_only = True
    for line in lines:
        if 'not pure fallback strategy anymore' in line:
            fallback_only = False
            break

    props['fallback_only'] = fallback_only

parser.add_function(parse_pure_fallback_strategy)

def parse_merging_for_symmetries_failures(content, props):
    merging_for_symmetries_attempts = props.get('merging_for_symmetries_attempts', 0)
    merging_for_symmetries_fail_shrinking = props.get('merging_for_symmetries_fail_shrinking', 0)
    merging_for_symmetries_fail_pruning = props.get('merging_for_symmetries_fail_pruning', 0)
    merging_for_symmetries_fail_any = props.get('merging_for_symmetries_fail_any', 0)

    merging_for_symmetries_success_ratio_shrinking = 1
    merging_for_symmetries_success_ratio_pruning = 1
    merging_for_symmetries_success_ratio_any = 1
    if merging_for_symmetries_attempts > 0:
        merging_for_symmetries_success_ratio_shrinking = 1 - float(merging_for_symmetries_fail_shrinking)/merging_for_symmetries_attempts
        merging_for_symmetries_success_ratio_pruning = 1 - float(merging_for_symmetries_fail_pruning)/merging_for_symmetries_attempts
        merging_for_symmetries_success_ratio_any = 1 - float(merging_for_symmetries_fail_any)/merging_for_symmetries_attempts

    props['merging_for_symmetries_success_ratio_shrinking'] = merging_for_symmetries_success_ratio_shrinking
    props['merging_for_symmetries_success_ratio_pruning'] = merging_for_symmetries_success_ratio_pruning
    props['merging_for_symmetries_success_ratio_any'] = merging_for_symmetries_success_ratio_any

parser.add_function(parse_merging_for_symmetries_failures)

parser.parse()

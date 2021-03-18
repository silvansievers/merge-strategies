#! /usr/bin/env python

from lab.parser import Parser
import re

parser = Parser()

parser.add_pattern('ms_avg_imperfect_shrinking', 'Average imperfect shrinking: (.+)', required=False, type=float)
parser.add_pattern('ms_not_exact_iteration', 'not perfect anymore in iteration (\d+)', required=False, type=int)
parser.add_pattern('ms_avg_pruning', 'Average relative pruning: (.+)', required=False, type=float)
parser.add_pattern('ms_tiebreaking_iterations', 'Iterations with merge tiebreaking: (\d+)', required=False, type=int)
parser.add_pattern('ms_tiebreaking_total', 'Total tiebreaking merge candidates: (\d+)', required=False, type=int)
parser.add_pattern('ms_max_int_abs_size', 'Maximum intermediate abstraction size: (\d+)', required=False, type=int)

def course_imperfect_shrinking(content, props):
    course = re.findall(r'Course of miss qualified states shrinking: \[(.*)\]', content)
    props['ms_course_imperfect_shrinking'] = course

parser.add_function(course_imperfect_shrinking)

def course_label_reduction(content, props):
    course = re.findall(r'Course of label reduction: \[(.*)\]', content)
    props['ms_course_label_reduction'] = course

parser.add_function(course_label_reduction)

def course_init_h_improvements(content, props):
    course = re.findall(r'Init h value improvements: \[(.*)\]', content)
    props['ms_init_h_improvements'] = course

parser.add_function(course_init_h_improvements)

def check_one_scc(content, props):
    one_scc = False
    for line in content.splitlines():
        if line == 'Only one single SCC':
            one_scc = True
            break
    props['ms_one_scc'] = one_scc

parser.add_function(check_one_scc)

def check_linear_order(content, props):
    linear_order = None
    for line in content.splitlines():
        if line == 'Linear merge order':
            linear_order = True
            break
        if line == 'Non-linear merge order':
            linear_order = False
            break
    props['ms_linear_order'] = linear_order

parser.add_function(check_linear_order)

def parse_merge_order(content, props):
    merge_order = re.findall(r'Merge order: \[(.*)\]', content)
    props['ms_merge_order'] = merge_order

parser.add_function(parse_merge_order)

def parse_course_pruning(content, props):
    course_pruning = re.findall(r'Relative pruning per iteration: \[(.*)\]', content)
    props['ms_course_pruning'] = course_pruning

parser.add_function(parse_course_pruning)

parser.parse()

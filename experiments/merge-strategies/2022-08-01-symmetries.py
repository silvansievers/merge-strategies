#! /usr/bin/env python3

import itertools
import os
from pathlib import Path
import subprocess

from lab.environments import LocalEnvironment, BaselSlurmEnvironment
from lab.reports import Attribute, arithmetic_mean, geometric_mean

from downward.reports.compare import ComparativeReport

import common_setup
from common_setup import IssueConfig, IssueExperiment

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ['DOWNWARD_BENCHMARKS']
REVISION='0f9cb8a9184a71c2c770538545b634dae1967328'
REVISIONS = [REVISION]
CONFIGS = [
    ## factored symmetries
    IssueConfig('b50k-cggl-largest-nonlinear-b60-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_symmetries(merge_tree=linear(variable_order=cg_goal_level,update_option=use_first),stabilize_transition_systems=false,max_bliss_iterations=infinity,bliss_call_time_limit=0,bliss_total_time_budget=60,symmetries_for_merging=largest,internal_merging=non_linear),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-dfp-largest-nonlinear-b60-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_symmetries(merge_selector=score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false,random_seed=2016)]),stabilize_transition_systems=false,max_bliss_iterations=infinity,bliss_call_time_limit=0,bliss_total_time_budget=60,symmetries_for_merging=largest,internal_merging=non_linear),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-miasm-dfp-largest-nonlinear-b60-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_symmetries(merge_tree=miasm(abstraction=miasm_merge_and_shrink,fallback_merge_selector=score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false,random_seed=2016)])),stabilize_transition_systems=false,max_bliss_iterations=infinity,bliss_call_time_limit=0,bliss_total_time_budget=60,symmetries_for_merging=largest,internal_merging=non_linear),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-rl-largest-nonlinear-b60-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_symmetries(merge_tree=linear(variable_order=reverse_level,update_option=use_first),stabilize_transition_systems=false,max_bliss_iterations=infinity,bliss_call_time_limit=0,bliss_total_time_budget=60,symmetries_for_merging=largest,internal_merging=non_linear),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-sbmiasm-largest-nonlinear-b60-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_symmetries(merge_selector=score_based_filtering(scoring_functions=[sf_miasm(shrink_strategy=shrink_bisimulation(greedy=false),max_states=50K,threshold_before_merge=1),total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false,random_seed=2016)]),stabilize_transition_systems=false,max_bliss_iterations=infinity,bliss_call_time_limit=0,bliss_total_time_budget=60,symmetries_for_merging=largest,internal_merging=non_linear),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
]

SUITE = common_setup.DEFAULT_OPTIMAL_SUITE
ENVIRONMENT = BaselSlurmEnvironment(
    email="silvan.sievers@unibas.ch",
    partition="infai_2",
    export=[],
    # paths obtained via:
    # module purge
    # module -q load Python/3.7.4-GCCcore-8.3.0
    # module -q load CMake/3.15.3-GCCcore-8.3.0
    # module -q load GCC/8.3.0
    # echo $PATH
    # echo $LD_LIBRARY_PATH
    setup='export PATH=/scicore/soft/apps/CMake/3.15.3-GCCcore-8.3.0/bin:/scicore/soft/apps/cURL/7.66.0-GCCcore-8.3.0/bin:/scicore/soft/apps/Python/3.7.4-GCCcore-8.3.0/bin:/scicore/soft/apps/XZ/5.2.4-GCCcore-8.3.0/bin:/scicore/soft/apps/SQLite/3.29.0-GCCcore-8.3.0/bin:/scicore/soft/apps/Tcl/8.6.9-GCCcore-8.3.0/bin:/scicore/soft/apps/ncurses/6.1-GCCcore-8.3.0/bin:/scicore/soft/apps/bzip2/1.0.8-GCCcore-8.3.0/bin:/scicore/soft/apps/binutils/2.32-GCCcore-8.3.0/bin:/scicore/soft/apps/GCCcore/8.3.0/bin:/infai/sieverss/repos/bin:/infai/sieverss/local:/export/soft/lua_lmod/centos7/lmod/lmod/libexec:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:$PATH\nexport LD_LIBRARY_PATH=/scicore/soft/apps/cURL/7.66.0-GCCcore-8.3.0/lib:/scicore/soft/apps/Python/3.7.4-GCCcore-8.3.0/lib:/scicore/soft/apps/libffi/3.2.1-GCCcore-8.3.0/lib64:/scicore/soft/apps/libffi/3.2.1-GCCcore-8.3.0/lib:/scicore/soft/apps/GMP/6.1.2-GCCcore-8.3.0/lib:/scicore/soft/apps/XZ/5.2.4-GCCcore-8.3.0/lib:/scicore/soft/apps/SQLite/3.29.0-GCCcore-8.3.0/lib:/scicore/soft/apps/Tcl/8.6.9-GCCcore-8.3.0/lib:/scicore/soft/apps/libreadline/8.0-GCCcore-8.3.0/lib:/scicore/soft/apps/ncurses/6.1-GCCcore-8.3.0/lib:/scicore/soft/apps/bzip2/1.0.8-GCCcore-8.3.0/lib:/scicore/soft/apps/binutils/2.32-GCCcore-8.3.0/lib:/scicore/soft/apps/zlib/1.2.11-GCCcore-8.3.0/lib:/scicore/soft/apps/GCCcore/8.3.0/lib64:/scicore/soft/apps/GCCcore/8.3.0/lib')

if common_setup.is_test_run():
    SUITE = IssueExperiment.DEFAULT_TEST_SUITE
    ENVIRONMENT = LocalEnvironment(processes=4)

exp = IssueExperiment(
    revisions=REVISIONS,
    configs=CONFIGS,
    environment=ENVIRONMENT,
)
exp.add_suite(BENCHMARKS_DIR, SUITE)

exp.add_parser(exp.EXITCODE_PARSER)
exp.add_parser(exp.TRANSLATOR_PARSER)
exp.add_parser(exp.SINGLE_SEARCH_PARSER)
exp.add_parser(exp.PLANNER_PARSER)

exp.add_parser('ms-parser.py')
exp.add_parser('merge-strategies-parser.py')
exp.add_parser('symmetries-parser.py')

exp.add_step('build', exp.build)
exp.add_step('start', exp.start_runs)
exp.add_fetcher(name='fetch')

extra_attributes=[
    Attribute('search_out_of_memory', absolute=True, min_wins=True),
    Attribute('search_out_of_time', absolute=True, min_wins=True),
    Attribute('ms_construction_time', absolute=False, min_wins=True, function=geometric_mean),
    Attribute('score_ms_construction_time', min_wins=False, digits=4),
    Attribute('ms_atomic_construction_time', absolute=False, min_wins=True, function=geometric_mean),
    Attribute('ms_abstraction_constructed', absolute=True, min_wins=False),
    Attribute('ms_atomic_fts_constructed', absolute=True, min_wins=False),
    Attribute('ms_out_of_memory', absolute=True, min_wins=True),
    Attribute('ms_out_of_time', absolute=True, min_wins=True),
    Attribute('ms_memory_delta', absolute=False, min_wins=True),
    Attribute('ms_reached_time_limit', absolute=False, min_wins=True),

    Attribute('ms_avg_imperfect_shrinking', absolute=False, min_wins=True, function=arithmetic_mean),
    Attribute('ms_course_imperfect_shrinking', absolute=True),
    Attribute('ms_course_label_reduction', absolute=True),
    Attribute('ms_init_h_improvements', absolute=False, min_wins=False),
    Attribute('ms_not_exact_iteration', absolute=False, min_wins=False),
    Attribute('ms_one_scc', absolute=True, min_wins=False),
    Attribute('ms_linear_order', absolute=True, min_wins=True),
    Attribute('ms_merge_order', absolute=True),
    Attribute('ms_course_pruning', absolute=True),
    Attribute('ms_avg_pruning', absolute=False, min_wins=False, function=arithmetic_mean),
    Attribute('ms_tiebreaking_iterations', absolute=True, min_wins=True),
    Attribute('ms_tiebreaking_total', absolute=True, min_wins=True),
    Attribute('ms_max_int_abs_size', absolute=False, min_wins=True, function=arithmetic_mean),

    Attribute('number_of_applied_symmetries', absolute=False, min_wins=False),
    Attribute('bliss_time_average', absolute=False, min_wins=True, function=arithmetic_mean),
    Attribute('bliss_time_median', absolute=False, min_wins=True, function=arithmetic_mean),
    Attribute('bliss_memory_out', absolute=True, min_wins=True),
    Attribute('bliss_timeout', absolute=True, min_wins=True),
    Attribute('bliss_total_calls', absolute=False, min_wins=False),
    Attribute('fallback_only', absolute=True, min_wins=True),
    Attribute('merging_for_symmetries_attempts', absolute=True, min_wins=False),
    Attribute('merging_for_symmetries_fail_shrinking', absolute=True, min_wins=True),
    Attribute('merging_for_symmetries_fail_pruning', absolute=True, min_wins=True),
    Attribute('merging_for_symmetries_fail_any', absolute=True, min_wins=True),
    Attribute('merging_for_symmetries_success_ratio_shrinking', absolute=False, min_wins=False, function=arithmetic_mean),
    Attribute('merging_for_symmetries_success_ratio_pruning', absolute=False, min_wins=False, function=arithmetic_mean),
    Attribute('merging_for_symmetries_success_ratio_any', absolute=False, min_wins=False, function=arithmetic_mean),
]
attributes = list(exp.DEFAULT_TABLE_ATTRIBUTES)
attributes.extend(extra_attributes)

exp.add_absolute_report_step(
    filter_algorithm=[
        '{}-b50k-cggl-largest-nonlinear-b60-time900'.format(REVISION),
        '{}-b50k-dfp-largest-nonlinear-b60-time900'.format(REVISION),
        '{}-b50k-miasm-dfp-largest-nonlinear-b60-time900'.format(REVISION),
        '{}-b50k-rl-largest-nonlinear-b60-time900'.format(REVISION),
        '{}-b50k-sbmiasm-largest-nonlinear-b60-time900'.format(REVISION),
    ],
    attributes=attributes,
)

OLD_REVISION='120ac87a77e6585786fe56fe1532ac82b8215007'
exp.add_fetcher(
    'data/2021-03-18-symmetries-eval',
    filter_algorithm=[
        '{}-b50k-cggl-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
        '{}-b50k-dfp-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
        '{}-b50k-miasm-dfp-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
        '{}-b50k-rl-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
        '{}-b50k-sbmiasm-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
    ],
    merge=True
)

report_name=f'{exp.name}-compare-old'
report_file=Path(exp.eval_dir) / f'{report_name}.html'
exp.add_report(
    ComparativeReport(
        attributes=attributes,
        algorithm_pairs=[
            ('{}-b50k-cggl-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
            '{}-b50k-cggl-largest-nonlinear-b60-time900'.format(REVISION)),
            ('{}-b50k-dfp-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
            '{}-b50k-dfp-largest-nonlinear-b60-time900'.format(REVISION)),
            ('{}-b50k-miasm-dfp-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
            '{}-b50k-miasm-dfp-largest-nonlinear-b60-time900'.format(REVISION)),
            ('{}-b50k-rl-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
            '{}-b50k-rl-largest-nonlinear-b60-time900'.format(REVISION)),
            ('{}-b50k-sbmiasm-largest-nonlinear-b60-partialmax-time900'.format(OLD_REVISION),
            '{}-b50k-sbmiasm-largest-nonlinear-b60-time900'.format(REVISION)),
        ],
    ),
    name=report_name,
    outfile=report_file,
)
exp.add_step('publish-comparison-report', subprocess.call, ['publish', report_file])

MIDDLE_REV='82a246792974adf880e4834762deef983d4389b1'
exp.add_fetcher(
    'data/2022-04-05-symmetries-eval',
    filter_algorithm=[
        '{}-b50k-cggl-largest-nonlinear-b60-time900'.format(MIDDLE_REV),
        '{}-b50k-dfp-largest-nonlinear-b60-time900'.format(MIDDLE_REV),
        '{}-b50k-miasm-dfp-largest-nonlinear-b60-time900'.format(MIDDLE_REV),
        '{}-b50k-rl-largest-nonlinear-b60-time900'.format(MIDDLE_REV),
        '{}-b50k-sbmiasm-largest-nonlinear-b60-time900'.format(MIDDLE_REV),
    ],
    merge=True
)

exp.add_comparison_table_step(revisions=[MIDDLE_REV, REVISION], attributes=attributes, name='compare-previous')

exp.run_steps()

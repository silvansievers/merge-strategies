#! /usr/bin/env python3

import itertools
import os

from lab.environments import LocalEnvironment, BaselSlurmEnvironment
from lab.reports import Attribute, arithmetic_mean, geometric_mean

from downward.reports.compare import ComparativeReport

import common_setup
from common_setup import IssueConfig, IssueExperiment

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ['DOWNWARD_BENCHMARKS']
REVISIONS = ['6b17a29694dbe2bc097a2b630931fe5e4c77733b']
CONFIGS = [
    ## scoring function collection
    IssueConfig('b50k-sbcollection-partialmax-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_stateless(merge_selector=score_based_weighted_sum(scoring_functions=[causal_connection,boolean_causal_connection,non_additivity,transitions_states_quotient(prefer_high=false),init_h(choice=improvement,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),max_fgh(fgh=f,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),avg_h(choice=improvement,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),goal_relevance_fine,num_variables,perfect_shrinking(shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),num_transitions,lr_opportunities,more_lr_opportunities,mutexes,total_order(atomic_ts_order=reverse_level,product_ts_order=old_to_new,atomic_before_product=true,random_seed=2016)])),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,partial_mas_method=maximum,main_loop_max_time=900),verbosity=silent)']),
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
]
attributes = list(exp.DEFAULT_TABLE_ATTRIBUTES)
attributes.extend(extra_attributes)

REVISION=REVISIONS[0]
exp.add_absolute_report_step(
    filter_algorithm=[
        '{}-b50k-sbcollection-partialmax-time900'.format(REVISION),
    ],
    attributes=attributes,
)

exp.run_steps()

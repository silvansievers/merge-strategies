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
OLD_REVISION='0f9cb8a9184a71c2c770538545b634dae1967328'
REVISION='1bca2102a7563b2fe93a3bb84a61f41394d20b52'
REVISIONS = [REVISION]
CONFIGS = [
    ## standard configs
    IssueConfig('b50k-cggl-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_precomputed(merge_tree=linear(variable_order=cg_goal_level)),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-dfp-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_stateless(merge_selector=score_based_filtering(scoring_functions=[goal_relevance(),dfp(),total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false)])),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-miasm-dfp-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_precomputed(merge_tree=miasm(abstraction=miasm_merge_and_shrink(),fallback_merge_selector=score_based_filtering(scoring_functions=[goal_relevance(),dfp(),total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false)]))),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-rl-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_precomputed(merge_tree=linear(variable_order=reverse_level)),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-sbmiasm-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_stateless(merge_selector=score_based_filtering(scoring_functions=[sf_miasm(shrink_strategy=shrink_bisimulation(greedy=false),max_states=50K,threshold_before_merge=1),total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false)])),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
    IssueConfig('b50k-sccdfp-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),merge_strategy=merge_sccs(order_of_sccs=topological,merge_selector=score_based_filtering(scoring_functions=[goal_relevance(),dfp(),total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false)])),max_states=50K,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),

    ## scoring function collection
    IssueConfig('b50k-sbcollection-time900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_stateless(merge_selector=score_based_weighted_sum(scoring_functions=[causal_connection(),boolean_causal_connection(),non_additivity(),transitions_states_quotient(prefer_high=false),init_h(choice=improvement,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),max_fgh(fgh=f,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),avg_h(choice=improvement,shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),goal_relevance_fine(),num_variables(),perfect_shrinking(shrink_strategy=shrink_bisimulation(greedy=false),max_states=50000,threshold_before_merge=1),num_transitions(),lr_opportunities(),more_lr_opportunities(),mutexes(),total_order(atomic_ts_order=reverse_level,product_ts_order=old_to_new,atomic_before_product=true,random_seed=2016)])),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,main_loop_max_time=900),verbosity=silent)']),
]

SUITE = common_setup.DEFAULT_OPTIMAL_SUITE
ENVIRONMENT = BaselSlurmEnvironment(
    partition="infai_3",
    email="silvan.sievers@unibas.ch",
    memory_per_cpu="3940M",
    cpus_per_task=2,
    export=["PATH"],
    # paths obtained via:
    # module purge
    # module -q load Python/3.10.4-GCCcore-11.3.0
    # module -q load GCC/11.3.0
    # module -q load CMake/3.23.1-GCCcore-11.3.0
    # echo $PATH
    # echo $LD_LIBRARY_PATH
    setup='export PATH=/scicore/soft/apps/CMake/3.23.1-GCCcore-11.3.0/bin:/scicore/soft/apps/libarchive/3.6.1-GCCcore-11.3.0/bin:/scicore/soft/apps/cURL/7.83.0-GCCcore-11.3.0/bin:/scicore/soft/apps/Python/3.10.4-GCCcore-11.3.0/bin:/scicore/soft/apps/OpenSSL/1.1/bin:/scicore/soft/apps/XZ/5.2.5-GCCcore-11.3.0/bin:/scicore/soft/apps/SQLite/3.38.3-GCCcore-11.3.0/bin:/scicore/soft/apps/Tcl/8.6.12-GCCcore-11.3.0/bin:/scicore/soft/apps/ncurses/6.3-GCCcore-11.3.0/bin:/scicore/soft/apps/bzip2/1.0.8-GCCcore-11.3.0/bin:/scicore/soft/apps/binutils/2.38-GCCcore-11.3.0/bin:/scicore/soft/apps/GCCcore/11.3.0/bin:/infai/sieverss/repos/bin:/infai/sieverss/local:/export/soft/lua_lmod/centos7/lmod/lmod/libexec:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:$PATH\nexport LD_LIBRARY_PATH=/scicore/soft/apps/libarchive/3.6.1-GCCcore-11.3.0/lib:/scicore/soft/apps/cURL/7.83.0-GCCcore-11.3.0/lib:/scicore/soft/apps/Python/3.10.4-GCCcore-11.3.0/lib:/scicore/soft/apps/OpenSSL/1.1/lib:/scicore/soft/apps/libffi/3.4.2-GCCcore-11.3.0/lib64:/scicore/soft/apps/GMP/6.2.1-GCCcore-11.3.0/lib:/scicore/soft/apps/XZ/5.2.5-GCCcore-11.3.0/lib:/scicore/soft/apps/SQLite/3.38.3-GCCcore-11.3.0/lib:/scicore/soft/apps/Tcl/8.6.12-GCCcore-11.3.0/lib:/scicore/soft/apps/libreadline/8.1.2-GCCcore-11.3.0/lib:/scicore/soft/apps/ncurses/6.3-GCCcore-11.3.0/lib:/scicore/soft/apps/bzip2/1.0.8-GCCcore-11.3.0/lib:/scicore/soft/apps/binutils/2.38-GCCcore-11.3.0/lib:/scicore/soft/apps/zlib/1.2.12-GCCcore-11.3.0/lib:/scicore/soft/apps/GCCcore/11.3.0/lib64')

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

exp.add_absolute_report_step(
    filter_algorithm=[f"{REVISION}-{config.nick}" for config in CONFIGS],
    attributes=attributes,
)

exp.add_fetcher(
    'data/2022-08-01-merge-strategies-eval',
    filter_algorithm=[f"{OLD_REVISION}-{config.nick}" for config in CONFIGS],
    merge=True
)

report_name=f'{exp.name}-compare-old'
report_file=Path(exp.eval_dir) / f'{report_name}.html'
exp.add_report(
    ComparativeReport(
        attributes=attributes,
        algorithm_pairs=[
            (f"{OLD_REVISION}-{config.nick}", f"{REVISION}-{config.nick}") for config in CONFIGS
        ],
    ),
    name=report_name,
    outfile=report_file,
)
exp.add_step('publish-comparison-report', subprocess.call, ['publish', report_file])

exp.run_steps()
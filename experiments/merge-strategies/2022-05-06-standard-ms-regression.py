#! /usr/bin/env python3

import itertools
import os
from pathlib import Path

from lab.environments import LocalEnvironment, BaselSlurmEnvironment
from lab.reports import Attribute, arithmetic_mean, geometric_mean

from downward.reports.compare import ComparativeReport

import common_setup
from common_setup import IssueConfig, IssueExperiment

DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.splitext(os.path.basename(__file__))[0]
BENCHMARKS_DIR = os.environ['DOWNWARD_BENCHMARKS']
REVISIONS = [
    'e966ed5db5c644402c28fe1ec54a7b39403b6b75', # FD main base before latest merge from main
    '5182b5afb1604b88af2270011affd6966f7a590b', # after issue1018 (PDB)
    '5430e040442147efebec24e7e404feff248d973b', # after issue995 (landmarks)
    '7516b8d7359049e284bf43fede87cfe1a33153dc', # after issue1007 (CEGAR PDBs)
    '40cb509cf40931db2b038558ce2a85e27bfda381', # after issue1008 (pattern generators)
    'bdfedac371f1daf3b80d9a74074c4172cabffb52', # after issue999 (landmarks)
    '8eba9ac9629f038086e038816fde745eb3a82cf9', # after issue1032 (RNG in VariableOrderFinder)
    'c39d79b53fc236743ff971b1ad42270ae15237e8', # after issue964 (logging)
    '63907c8b0291d8eca910bf6ce5a5e2075b74f261', # after issue467 (landmarks)
    'f277c0c482511a8f877e57d973dc6139cc1f17e7', # after issue1041 (landmarks)
    'ec96fbdb53decee03658761e2945f33cf50b0f21', # after issue998 (landmarks)
    '75a889eadc22aa8fea97c406364f0f0c64131671', # after issue983 (delete relaxation)
    '6769e68c6a728ab2a76cc12e2e6ae066f7fd289e', # after issue937 (landmarks)
    '2a9cc5a5a6870031808c4e896393c1c82d9289b2', # after issue1042 (limited pruning)
    'b772e1c3454c4a750498b943fcb553a3e3634c37', # after issue1043 (logging) = FD 21.12
]
CONFIGS = [
    IssueConfig('dfp-b50k-t900', ['--search', 'astar(merge_and_shrink(merge_strategy=merge_stateless(merge_selector=score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order])),shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,main_loop_max_time=900))']),
    # IssueConfig('rl-b50k-t900', ['--search', 'astar(merge_and_shrink(merge_strategy=merge_precomputed(merge_tree=linear(variable_order=reverse_level)),shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,main_loop_max_time=900))']),
    # IssueConfig('sccs-dfp-b50k-t900', ['--search', 'astar(merge_and_shrink(shrink_strategy=shrink_bisimulation(greedy=false),merge_strategy=merge_sccs(order_of_sccs=topological,merge_selector=score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order(atomic_ts_order=reverse_level,product_ts_order=new_to_old,atomic_before_product=false)])),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,main_loop_max_time=900))']),
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
]
attributes = list(exp.DEFAULT_TABLE_ATTRIBUTES)
attributes.extend(extra_attributes)

exp.add_absolute_report_step(attributes=attributes)

for i in range(1, len(REVISIONS)):
    rev1 = REVISIONS[i-1]
    rev2 = REVISIONS[i]
    report_name=f'{exp.name}-compare-{rev1}-{rev2}'
    report_file=Path(exp.eval_dir) / f'{report_name}.html'
    exp.add_report(
        ComparativeReport(
            attributes=attributes,
            algorithm_pairs=[
                ('{}-dfp-b50k-t900'.format(rev1),
                '{}-dfp-b50k-t900'.format(rev2)),
            ],
        ),
        name=report_name,
        outfile=report_file,
    )

exp.run_steps()

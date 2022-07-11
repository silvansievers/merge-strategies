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
    '6769e68c6a728ab2a76cc12e2e6ae066f7fd289e', # after issue937 (landmarks)
    '2a9cc5a5a6870031808c4e896393c1c82d9289b2', # after issue1042 (limited pruning)
]
CONFIGS = [
    IssueConfig('dfp-b50k-t900', ['--search', 'astar(merge_and_shrink(merge_strategy=merge_stateless(merge_selector=score_based_filtering(scoring_functions=[goal_relevance,dfp,total_order])),shrink_strategy=shrink_bisimulation(greedy=false),label_reduction=exact(before_shrinking=true,before_merging=false),max_states=50000,threshold_before_merge=1,main_loop_max_time=900))']),
    IssueConfig('lmcut', ['--search', 'astar(lmcut())']),
    IssueConfig('cpdbs', ['--search', 'astar(cpdbs(multiple_cegar(max_pdb_size=1000000,max_collection_size=10000000,pattern_generation_max_time=infinity,total_max_time=3,stagnation_limit=20,blacklist_trigger_percentage=0.75,enable_blacklist_on_stagnation=true,random_seed=2018,verbosity=normal,use_wildcard_plans=false)),verbosity=normal)']),
    IssueConfig('bjolp', ['--evaluator', 'lmc=lmcount(lm_merged([lm_rhw(),lm_hm(m=1)]),admissible=true)', '--search', 'astar(lmc,lazy_evaluator=lmc)']),
    IssueConfig('blind', ['--search', 'astar(blind())']),
]

SUITE = common_setup.DEFAULT_OPTIMAL_SUITE
ENVIRONMENT = BaselSlurmEnvironment(
    email="silvan.sievers@unibas.ch",
    partition="infai_2",
    export=[],
    # paths obtained via:
    # module purge
    # module -q load CMake/3.20.1-GCCcore-10.3.0
    # module -q load Clang/12.0.1-GCCcore-10.3.0
    # echo $PATH
    # echo $LD_LIBRARY_PATH
    # use export CC=clang and export CXX=clang++ before compiling!
    setup='export PATH=/scicore/soft/apps/Clang/12.0.1-GCCcore-10.3.0/bin:/scicore/soft/apps/Z3/4.8.11-GCCcore-10.3.0/bin:/scicore/soft/apps/hwloc/2.4.1-GCCcore-10.3.0/sbin:/scicore/soft/apps/hwloc/2.4.1-GCCcore-10.3.0/bin:/scicore/soft/apps/libxml2/2.9.10-GCCcore-10.3.0/bin:/scicore/soft/apps/numactl/2.0.14-GCCcore-10.3.0/bin:/scicore/soft/apps/binutils/2.36.1-GCCcore-10.3.0/bin:/scicore/soft/apps/CMake/3.20.1-GCCcore-10.3.0/bin:/scicore/soft/apps/libarchive/3.5.1-GCCcore-10.3.0/bin:/scicore/soft/apps/XZ/5.2.5-GCCcore-10.3.0/bin:/scicore/soft/apps/cURL/7.76.0-GCCcore-10.3.0/bin:/scicore/soft/apps/bzip2/1.0.8-GCCcore-10.3.0/bin:/scicore/soft/apps/ncurses/6.2-GCCcore-10.3.0/bin:/scicore/soft/apps/GCCcore/10.3.0/bin:/infai/sieverss/repos/bin:/infai/sieverss/local:/export/soft/lua_lmod/centos7/lmod/lmod/libexec:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin\nexport LD_LIBRARY_PATH=/scicore/soft/apps/Clang/12.0.1-GCCcore-10.3.0/lib:/scicore/soft/apps/Z3/4.8.11-GCCcore-10.3.0/lib:/scicore/soft/apps/GMP/6.2.1-GCCcore-10.3.0/lib:/scicore/soft/apps/hwloc/2.4.1-GCCcore-10.3.0/lib:/scicore/soft/apps/libpciaccess/0.16-GCCcore-10.3.0/lib:/scicore/soft/apps/libxml2/2.9.10-GCCcore-10.3.0/lib:/scicore/soft/apps/numactl/2.0.14-GCCcore-10.3.0/lib:/scicore/soft/apps/binutils/2.36.1-GCCcore-10.3.0/lib:/scicore/soft/apps/libarchive/3.5.1-GCCcore-10.3.0/lib:/scicore/soft/apps/XZ/5.2.5-GCCcore-10.3.0/lib:/scicore/soft/apps/cURL/7.76.0-GCCcore-10.3.0/lib:/scicore/soft/apps/bzip2/1.0.8-GCCcore-10.3.0/lib:/scicore/soft/apps/zlib/1.2.11-GCCcore-10.3.0/lib:/scicore/soft/apps/ncurses/6.2-GCCcore-10.3.0/lib:/scicore/soft/apps/GCCcore/10.3.0/lib64')

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

exp.add_comparison_table_step(attributes=attributes)

exp.run_steps()

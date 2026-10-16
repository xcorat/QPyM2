"""
Grid parser utility for toy batch configurations.

Parses model_grid from toy config YAML files and outputs all combinations
of grid parameters for parallel batch job submission.

Usage:
    python -m qpym2.utils.grid_parser --config batch
    python -m qpym2.utils.grid_parser --config batch --format json
"""

import argparse
import json
import sys
from typing import List, Dict, Any
from qpym2.cfg.toy_config import load_toy_config


def parse_grid(config_name: str = 'batch') -> List[Dict[str, Any]]:
    """
    Parse model grid from toy config and return all combinations.
    
    Args:
        config_name: Name of toy config to load (default: 'batch')
        
    Returns:
        List of dicts, each containing one grid point:
        {'co_cut': str, 'ndbd_cutoff': str, 'co_shift': float}
    
    Examples:
        >>> combos = parse_grid('batch')
        >>> len(combos)  # 3 co_cuts × 3 ndbd_cutoff × 7 co_shift
        63
        >>> combos[0]
        {'co_cut': 'gamma10_v150', 'ndbd_cutoff': 'u1740', 'co_shift': 0.0}
    """
    cfg = load_toy_config(config_name)
    
    if 'model_grid' not in cfg:
        raise ValueError(
            f"Toy config '{config_name}' has no 'model_grid' section.\n"
            f"Expected: model_grid: {{co_cuts: [...], ndbd_cutoff: [...], co_shift: [...]}}"
        )
    
    grid = cfg['model_grid']
    
    # Validate required fields
    required = ['co_cuts', 'ndbd_cutoff', 'co_shift']
    missing = [field for field in required if field not in grid]
    if missing:
        raise ValueError(
            f"model_grid missing required fields: {missing}\n"
            f"Required: {required}"
        )
    
    # Build Cartesian product
    combinations = []
    for co_cut in grid['co_cuts']:
        for ndbd_cutoff in grid['ndbd_cutoff']:
            for co_shift in grid['co_shift']:
                combinations.append({
                    'co_cut': co_cut,
                    'ndbd_cutoff': ndbd_cutoff,
                    'co_shift': co_shift
                })
    
    return combinations


def format_csv(combinations: List[Dict[str, Any]]) -> str:
    """Format grid combinations as CSV (one per line)."""
    lines = []
    for combo in combinations:
        # Handle None/null values
        co_cut = combo['co_cut'] if combo['co_cut'] is not None else 'null'
        ndbd_cutoff = combo['ndbd_cutoff'] if combo['ndbd_cutoff'] is not None else 'null'
        co_shift = combo['co_shift'] if combo['co_shift'] is not None else 'null'
        lines.append(f"{co_cut},{ndbd_cutoff},{co_shift}")
    return '\n'.join(lines)


def format_json(combinations: List[Dict[str, Any]]) -> str:
    """Format grid combinations as JSON array."""
    return json.dumps(combinations, indent=2)


def main():
    """CLI entry point for grid parser."""
    parser = argparse.ArgumentParser(
        description='Parse toy config model_grid and output all combinations.'
    )
    parser.add_argument(
        '--config', '-c',
        default='batch',
        help='Toy config name (default: batch)'
    )
    parser.add_argument(
        '--format', '-f',
        choices=['csv', 'json'],
        default='csv',
        help='Output format (default: csv)'
    )
    parser.add_argument(
        '--count',
        action='store_true',
        help='Only print the count of combinations'
    )
    
    args = parser.parse_args()
    
    try:
        combinations = parse_grid(args.config)
        
        if args.count:
            print(len(combinations))
        elif args.format == 'csv':
            print(format_csv(combinations))
        elif args.format == 'json':
            print(format_json(combinations))
        
        return 0
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

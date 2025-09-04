#!/usr/bin/env python3
"""
RNAcentral Top Hit Extractor

A Click command-line tool to extract the top hit from RNAcentral JSON search results.
"""

import json
import sys
from pathlib import Path
import click
from typing import Dict, Any, Optional

from pygments.lexer import default


@click.command()
@click.argument('json_file', type=click.Path(exists=True, path_type=Path))
@click.option('--output', '-o', type=click.Path(path_type=Path),
              help='Output file path (default: stdout)')
@click.option('--format', 'output_format', type=click.Choice(['json', 'tsv', 'summary']),
              default='json', help='Output format')
@click.option('--fields', '-f', multiple=True,
              help='Specific fields to extract (can be used multiple times)')
@click.option('--species', '-s', multiple=False, default="any",
              help='Extract top hit from a specific species')
@click.option('--quiet', '-q', is_flag=True, help='Suppress informational messages')
def extract_top_hit(json_file: Path, output: Optional[Path], output_format: str,
                    fields: tuple, species:str, quiet: bool):
    """
    Extract the top hit from RNAcentral JSON search results.

    JSON_FILE: Path to the RNAcentral JSON results file

    Examples:

        # Extract top hit as JSON
        rnacentral-extract results.json

        # Extract specific fields as TSV
        rnacentral-extract results.json --format tsv --fields rnacentral_id --fields rna_type --fields e_value

        # Save to file with summary format
        rnacentral-extract results.json -o top_hit.txt --format summary
    """

    try:
        # Load JSON data
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Check if results exist
        if 'results' not in data:
            click.echo(f"Error: No 'results' field found in {json_file}", err=True)
            sys.exit(1)

        if not data['results']:
            click.echo(f"Error: No results found in {json_file}", err=True)
            sys.exit(1)

        # Get top hit
        if species == 'any':
            top_hit = data['results'][0]
        else:
            for hit in data['results']:
                if species in hit['description']:
                    top_hit = hit
                    break
        # Process output based on format
        output_text = format_output(top_hit, output_format, fields)

        # Write output
        if output:
            with open(output, 'w') as f:
                f.write(output_text)
            if not quiet:
                click.echo(f"Top hit extracted to: {output}")
        else:
            click.echo(output_text)

        # Print summary info to stderr if not quiet
        if not quiet and output:
            total_results = len(data['results'])
            click.echo(f"Processed {total_results} results from {json_file}", err=True)

    except FileNotFoundError:
        click.echo(f"Error: File {json_file} not found", err=True)
        sys.exit(1)
    except json.JSONDecodeError as e:
        click.echo(f"Error: Invalid JSON in {json_file}: {e}", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def format_output(hit: Dict[Any, Any], output_format: str, fields: tuple) -> str:
    """Format the top hit according to the specified output format."""

    if output_format == 'json':
        if fields:
            # Extract only specified fields
            filtered_hit = {}
            for field in fields:
                if field in hit:
                    filtered_hit[field] = hit[field]
                elif field in hit.get('fields', {}):
                    filtered_hit[field] = hit['fields'][field]
                else:
                    filtered_hit[field] = None
            return json.dumps(filtered_hit, indent=2)
        else:
            return json.dumps(hit, indent=2)

    elif output_format == 'tsv':
        if fields:
            values = []
            for field in fields:
                if field in hit:
                    value = hit[field]
                elif field in hit.get('fields', {}):
                    value = hit['fields'][field]
                    # Handle lists in fields
                    if isinstance(value, list):
                        value = ';'.join(map(str, value)) if value else ''
                else:
                    value = ''
                values.append(str(value))

            # Header and data
            header = '\t'.join(fields)
            data = '\t'.join(values)
            return f"{header}\n{data}"
        else:
            # Default fields for TSV
            default_fields = ['rnacentral_id', 'description', 'score', 'e_value', 'identity', 'query_coverage']
            return format_output(hit, 'tsv', default_fields)

    elif output_format == 'summary':
        rna_type = hit.get('fields', {}).get('rna_type', ['Unknown'])
        rna_type = rna_type[0] if isinstance(rna_type, list) and rna_type else 'Unknown'

        summary = f"""RNAcentral Top Hit Summary
==========================
RNAcentral ID: {hit.get('rnacentral_id', 'N/A')}
Description: {hit.get('description', 'N/A')}
RNA Type: {rna_type}
Score: {hit.get('score', 'N/A')}
E-value: {hit.get('e_value', 'N/A')}
Identity: {hit.get('identity', 'N/A'):.2f}% (if available)
Query Coverage: {hit.get('query_coverage', 'N/A'):.2f}% (if available)
Target Length: {hit.get('target_length', 'N/A')} nt
Alignment Length: {hit.get('alignment_length', 'N/A')} nt

Expert Databases: {'; '.join(hit.get('fields', {}).get('expert_db', [])) or 'N/A'}
Gene Names: {'; '.join(hit.get('fields', {}).get('gene', [])) or 'N/A'}
"""

        if hit.get('alignment'):
            summary += f"\nAlignment:\n{hit['alignment']}\n"

        return summary

    else:
        raise ValueError(f"Unknown output format: {output_format}")


@click.command()
@click.argument('directory', type=click.Path(exists=True, path_type=Path))
@click.option('--pattern', '-p', default='*.json', help='File pattern to match')
@click.option('--output-dir', '-d', type=click.Path(path_type=Path),
              help='Output directory (default: same as input)')
@click.option('--format', 'output_format', type=click.Choice(['json', 'tsv', 'summary']),
              default='tsv', help='Output format')
@click.option('--suffix', '-s', default='_top_hit', help='Suffix for output files')
def batch_extract(directory: Path, pattern: str, output_dir: Optional[Path],
                  output_format: str, suffix: str):
    """
    Batch extract top hits from multiple RNAcentral JSON files.

    DIRECTORY: Directory containing RNAcentral JSON result files
    """

    if output_dir is None:
        output_dir = directory
    else:
        output_dir.mkdir(exist_ok=True)

    json_files = list(directory.glob(pattern))

    if not json_files:
        click.echo(f"No files matching pattern '{pattern}' found in {directory}", err=True)
        sys.exit(1)

    click.echo(f"Processing {len(json_files)} files...")

    processed = 0
    failed = 0

    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)

            if 'results' not in data or not data['results']:
                click.echo(f"Skipping {json_file.name}: No results found", err=True)
                failed += 1
                continue

            top_hit = get_top_hit_by_species(data['results'], species)

            if not top_hit:
                species_msg = f" for species '{species}'" if species != 'any' else ""
                click.echo(f"Skipping {json_file.name}: No hits found{species_msg}", err=True)
                failed += 1
                continue

            # Determine output file extension
            ext_map = {'json': '.json', 'tsv': '.tsv', 'summary': '.txt'}
            output_file = output_dir / f"{json_file.stem}{suffix}{ext_map[output_format]}"

            output_text = format_output(top_hit, output_format, ())

            with open(output_file, 'w') as f:
                f.write(output_text)

            processed += 1

        except Exception as e:
            click.echo(f"Error processing {json_file.name}: {e}", err=True)
            failed += 1

    click.echo(f"Completed: {processed} processed, {failed} failed")


@click.group()
@click.version_option(version='1.0.0')
def cli():
    """RNAcentral JSON Results Extractor"""
    pass


# Add commands to the group
cli.add_command(extract_top_hit)
cli.add_command(batch_extract)

if __name__ == '__main__':
    cli()
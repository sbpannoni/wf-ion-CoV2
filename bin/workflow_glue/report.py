#!/usr/bin/env python
"""Create report file."""

import json

from aplanat import bars, lines, report
from aplanat.components import bcfstats, nextclade
from aplanat.components import simple as scomponents
from aplanat.util import Colors
from bokeh.layouts import gridplot, layout
from bokeh.models import Panel, Range1d, Tabs
import pandas as pd
from pandas.api import types as pd_types
import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101

# Define categorical types
CATEGORICAL = pd_types.CategoricalDtype(ordered=True)


def read_files(summaries, **kwargs):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, **kwargs))
    return pd.concat(dfs)


def output_json(df, consensus_fasta):
    """Read depth stats df and create JSON output."""
    grouped_by_sample = df.groupby('sample_name')
    all_json = {}
    for sample in grouped_by_sample.groups.keys():
        group_by_primer = grouped_by_sample.get_group(
            sample).groupby('primer_set')
        rg1 = group_by_primer.get_group(1).reset_index()
        rg2 = group_by_primer.get_group(2).reset_index()
        newdf = pd.DataFrame()
        newdf['start'] = rg1['pos'] - 10
        newdf['start'] = newdf['start'].clip(lower=0).astype(int)
        newdf['end'] = rg1['pos'] + 10
        newdf['end'] = newdf['end'].astype(int)
        newdf['fwd'] = rg1['depth_fwd'] + rg2['depth_fwd']
        newdf['rev'] = rg1['depth_rev'] + rg2['depth_rev']
        newdf['rg1'] = rg1['depth']
        newdf['rg2'] = rg2['depth']
        # preserve data type by converting column by column
        final = list(list(x) for x in zip(*(
            newdf[x].values.tolist() for x in newdf.columns)))
        all_json[sample] = final
    final_json = {'data': []}
    # parse the consensus fasta to get extra info required
    with pysam.FastxFile(consensus_fasta) as fh:
        for entry in fh:
            all_json[entry.name][-1][1] = len(entry.sequence)
            final_json['data'].append({
                'barcode': entry.name,
                'chromosome': entry.comment,
                'seqlen': len(entry.sequence),
                'ncount': (entry.sequence).count('N'),
                'coverage': all_json[entry.name]
            })
    return final_json


def main(args):
    """Run entry point."""
    logger = get_named_logger("report")

    report_doc = report.WFReport(
        "SARS-CoV-2 Ion Torrent Ampliseq Sequencing Report", "wf-ion-CoV2",
        revision=args.revision, commit=args.commit)

    # Get samples and types into something we can use
    # Not currently in use but here for future use
    # sample_types = [{'sample': sample, 'type': type}
    #       for sample, type in zip(args.samples, args.types)]

    section = report_doc.add_section()
    section.markdown('''
### Ion Torrent Ampliseq Analysis Status

The panel below lists samples which failed to produce results. Samples not
listed here were analysed successfully, but may still contain inconclusive or
invalid results. See the following sections for further information.
''')
    status = pd.read_csv(args.status, sep='\t')
    failed = status.loc[status['pass'] == 0]
    if len(failed) == 0:
        fail_list = "All samples analysed successfully"
    else:
        fail_list = failed['sample'].str.cat(sep=', ')
    section.markdown("""
```{}```
""".format(fail_list))
    fail_percentage = int(100 * len(failed) / len(status))
    classes = ['Success', 'Analysis Failed']
    values = [100 - fail_percentage, fail_percentage]
    colors = ['#54B8B1', '#EF4135']
    plot = bars.single_hbar(
        values, classes, colors,
        title="Completed analyses",
        x_axis_label="%age Samples")
    plot.x_range = Range1d(0, 140)
    section.plot(plot)

    if not args.hide_coverage:
        section = report_doc.add_section()
        section.markdown('''
### Genome coverage

Plots below indicate depth of coverage coloured by Ampliseq primer pool.
iVar variant calling requires a minimum depth of at least {}X; positions
below this threshold are masked as N in the consensus.

***NB: Depth axes are not tied between samples — do not compare absolute
values across plots.***
'''.format(args.report_depth))

        # depth summary by amplicon pool
        df = read_files(
            args.depths, sep="\t", converters={'sample_name': str})
        epi2me_json = output_json(df, args.consensus_fasta)
        json_object = json.dumps(epi2me_json, indent=4, separators=(',', ':'))
        json_file = open("artic.json", "a")
        json_file.write(json_object)
        json_file.close()
        plots_pool = list()
        plots_orient = list()
        plots_combined = list()
        depth_lim = args.report_depth
        for sample in sorted(df['sample_name'].unique()):
            bc = df['sample_name'] == sample
            depth = df[bc].groupby('pos').sum().reset_index()
            depth_thresh = 100*(
                depth['depth'] >= depth_lim).sum() / len(depth['depth'])
            depth_mean = depth['depth'].mean()

            # total depth plot
            # plot line just to get aplanat niceities
            p = lines.line(
                [depth['pos']], [depth['depth']], colors=[Colors.cerulean],
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth',
            )
            p.varea(
                x=depth['pos'], y1=0.1, y2=depth['depth'],
                fill_color=Colors.cerulean)
            plots_combined.append(p)

            # fwd/rev
            xs = [depth['pos'], depth['pos']]
            ys = [depth['depth_fwd'], depth['depth_rev']]
            names = ['fwd', 'rev']
            colors = [Colors.dark_gray, Colors.verdigris]

            p = lines.line(
                xs, ys, colors=colors, names=names,
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth')
            for x, y, name, color in zip(xs, ys, names, colors):
                p.varea(
                    x=x, y1=0, y2=y, legend_label=name,
                    fill_color=color, alpha=0.7,
                    muted_color=color, muted_alpha=0.2)
            p.legend.click_policy = 'mute'
            plots_orient.append(p)

            # primer set plot
            pset = df['primer_set']
            xs = [df.loc[(pset == i) & bc]['pos'] for i in (1, 2)]
            ys = [df.loc[(pset == i) & bc]['depth'] for i in (1, 2)]
            names = ['pool-1', 'pool-2']
            colors = [Colors.light_cornflower_blue, Colors.feldgrau]

            p = lines.line(
                xs, ys, colors=colors, names=names,
                title="{}: {:.0f}X, {:.1f}% > {}X".format(
                    sample, depth_mean, depth_thresh, depth_lim),
                height=250, width=400,
                x_axis_label='position', y_axis_label='depth')
            for x, y, name, color in zip(xs, ys, names, colors):
                p.varea(
                    x=x, y1=0, y2=y, legend_label=name,
                    fill_color=color, alpha=0.7,
                    muted_color=color, muted_alpha=0.2)
            p.legend.click_policy = 'mute'
            plots_pool.append(p)

        tab1 = Panel(
            child=gridplot(plots_combined, ncols=3), title="Coverage Plot")
        tab2 = Panel(
            child=gridplot(plots_pool, ncols=3), title="By amplicon pool")
        tab3 = Panel(
            child=gridplot(plots_orient, ncols=3), title="By read orientation")
        cover_panel = Tabs(tabs=[tab1, tab2, tab3])
        section.plot(cover_panel)

    # canned VCF stats report component
    if not args.hide_variants:
        section = report_doc.add_section()
        bcfstats.full_report(args.bcftools_stats, report=section)

    # NextClade analysis
    if args.nextclade is not None:
        section = report_doc.add_section(
            section=nextclade.NextClade(args.nextclade))
        section.markdown(
            "The table shows errors, warnings or failed genes per sample:")
        error_df = pd.read_csv(args.nextclade_errors, sep=";").fillna('None')
        error_df = error_df[
            ['index', 'seqName', 'clade', 'failedCdses', 'warnings', 'errors']]
        error_df.rename(
            columns={
                'seqName': 'Sample Name',
                'warnings': 'Warnings',
                'errors': 'Errors',
                'failedCdses': 'Failed CDS'}, inplace=True)
        section.table(error_df, index=False)

    # Pangolin analysis
    if args.pangolin is not None:
        section = report_doc.add_section()
        section.markdown('''
### Lineage

The table below reports the lineage of each sample as calculated by
[pangolin](https://github.com/cov-lineages/pangolin).

''')
        section.table(pd.read_csv(args.pangolin), index=False)

    # Genotyping
    if args.genotypes is not None:
        section = report_doc.add_section()
        section.markdown('''
### Genotyping

The table below lists whether candidate variants were determined to exist
within each sample.

Variants are pre-filtered by iVar using minimum quality (default 20) and
minimum allele frequency (default 0.25). Positions with coverage below the
minimum depth threshold are reported as `n/a`.
''')
        df = read_files(args.genotypes, sep=',')
        df = df[[
            'Sample', 'Result', 'Date Tested', 'Lab ID', 'testKit',
            'CH1-Target', 'CH1-Result', 'CH1-Conf']]
        df = df.sort_values(by=['Sample', 'CH1-Target'], ascending=True)
        section.table(df, index=False)

    # Versions and params
    section = report_doc.add_section(
        section=scomponents.version_table(args.versions))
    section = report_doc.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report_doc.write(args.output)
    logger.info(f"Written report to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")

    parser.add_argument(
        "status",
        help="artic status file")
    parser.add_argument(
        "output",
        help="Report output filename")
    parser.add_argument(
        "--pangolin",
        help="pangolin CSV output file")
    parser.add_argument(
        "--depths", nargs='+', required=True,
        help="Depth summary files")
    parser.add_argument(
        "--nextclade",
        help="nextclade json output file")
    parser.add_argument(
        "--nextclade_errors",
        help="nextclade error csv")
    parser.add_argument(
        "--bcftools_stats", nargs='+', required=True,
        help="Outputs from bcftools stats")
    parser.add_argument(
        "--genotypes", nargs='+', required=False,
        help="Genotyping summary files")
    parser.add_argument(
        "--report_depth", default=100, type=int,
        help=(
            "Depth at which to provide a coverage statistics, "
            "e.g. 76% of genome covered at `report_depth`"))
    parser.add_argument(
        "--hide_coverage", action="store_true",
        help="Do not display coverage plots in report.")
    parser.add_argument(
        "--hide_variants", action="store_true",
        help="Do not display variant summary in report.")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--consensus_fasta",
        help="Fasta of all conesusus seqeunces")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")

    return parser

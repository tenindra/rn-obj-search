#!/usr/bin/env python3
"""Draws the parameter test query results (runParameterExperiments) of each test as one combined figure, both as they are and
relative to each query type's default setting.

For each test, writes <network>_parameter_<test>_<query>_query_relative.txt next to the stats file of each query type of
COMBINED_QUERIES, in the same format, with the query time (column 7) divided by that of the query type's default setting, the
test defaults of the dataset config (par_<query>_defaults, e.g. fanout=8,numlmks=1,lmktype=border_minmax,lambda=256):
  density_vs_<p> (one block per value of p, rows by density): at each density, by the default value of p's row at that density,
                 so the default series is flat at 1;
  <p>_vs_<q>     (one block per value of q, rows by p):        all by one value, the row of the default p and q, so that the blocks
                 stay comparable.
A query type whose default setting is not among the tested values is reported and its test skipped. As each query type has its
own defaults, each panel has its own reference.
Each figure puts the query types side by side (gnuplot multiplot), with one key above them, as the series are the same in every
panel: <fig_dir>/<network>_parameter_<test>_query_combined.eps with the query times, each panel with the y axis of its own
figure, and <...>_query_relative_combined.eps with the relative times, on linear y axes labelled "Relative Query Time" (the
reference of each panel is printed for the caption). The panels reuse the gnuplot
scripts generateFigure saved next to the original data (<...>.gp, so runParameterExperiments must have drawn its figures), so
titles, colours and tics are the same. The original files are not modified. Every network with parameter test results in
stats_dir (globalVariables) is drawn.
"""

import os
import re
import subprocess
import sys

TESTS = ["density_vs_numlmks", "density_vs_lmktype", "density_vs_lambda", "density_vs_fanout",
         "fanout_vs_lambda", "numlmks_vs_lambda", "numlmks_vs_lmktype"]
# Panels of the combined figures, left to right, with their titles (gnuplot enhanced text: an italic k)
COMBINED_QUERIES = [("aknn", "A{/:Italic k}NN"), ("fkn", "{/:Italic k}FN"), ("range", "Range")]
PARAMETERS_COL = 4  # 1-based, as in the stats files and gnuplot: the key=value list of the tested setting
TIME_COL = 7
DENSITY_COL = 9
SYMBOLS = {"numlmks": "m", "fanout": "b", "lambda": "λ"}  # The parameters' names in the printed references

USAGE = """Usage: python3 generateCombinedParameterFigures.py
Draws the parameter test query results (runParameterExperiments) of each network with the AkNN, FkN and range panels of each
test in one figure, as they are (fig_dir/<network>_parameter_<test>_query_combined.eps) and relative to each query type's default
setting, par_<query>_defaults in the dataset config (fig_dir/<network>_parameter_<test>_query_relative_combined.eps)."""


def read_blocks(path):
    """Blocks of rows (lists of fields) separated by blank lines, as gnuplot's index reads them."""
    blocks, block = [], []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.strip() == "":
                if block:
                    blocks.append(block)
                    block = []
            else:
                block.append(line.split("\t"))
    if block:
        blocks.append(block)
    return blocks


def write_blocks(path, blocks):
    with open(path, "w") as f:
        f.write("\n\n\n".join("\n".join("\t".join(row) for row in block) for block in blocks) + "\n")


def same_value(a, b):
    try:
        return float(a) == float(b)
    except ValueError:
        return a == b


def has_defaults(row, defaults, parameters):
    """Whether the setting of a stats row has the default value of each of the parameters."""
    setting = dict(pair.split("=", 1) for pair in row[PARAMETERS_COL - 1].split(",") if "=" in pair)
    return all(p in setting and same_value(setting[p], defaults[p]) for p in parameters)


def reference_parameters(test):
    """The parameters whose default values select the reference rows: the series parameter of a density test (one row per
    density), the x and series parameters otherwise (one row)."""
    x_param, series_param = test.split("_vs_")
    return [series_param] if x_param == "density" else [x_param, series_param]


def relative_blocks(blocks, test, defaults):
    """The blocks with each query time divided by the default setting's, or None if the default setting was not tested."""
    ref_rows = [row for block in blocks for row in block if has_defaults(row, defaults, reference_parameters(test))]
    if test.startswith("density_vs_"):
        ref = {row[DENSITY_COL - 1]: float(row[TIME_COL - 1]) for row in ref_rows}
        if any(row[DENSITY_COL - 1] not in ref for block in blocks for row in block):
            return None
        value = lambda row: ref[row[DENSITY_COL - 1]]
    else:
        # One reference for the whole figure: dividing each block by its own default row would hide the differences between blocks
        if len(ref_rows) != 1:
            return None
        value = lambda row: float(ref_rows[0][TIME_COL - 1])
    return [[row[:TIME_COL - 1] + ["%.4f" % (float(row[TIME_COL - 1]) / value(row))] + row[TIME_COL:] for row in block]
            for block in blocks]


def reference_label(test, defaults):
    """The default setting a panel is relative to, as the parameters' values (e.g. b=8, or λ=256, b=4)."""
    return ", ".join(defaults[p] if p == "lmktype" else "%s=%s" % (SYMBOLS[p], defaults[p])
                     for p in reversed(reference_parameters(test)))


def combined_figure(stats_dir, fig_dir, network, test, relative, defaults):
    """Draws the test's results of each query type of COMBINED_QUERIES as one figure with one key above the panels, as they are or
    (relative) converted relative to each query type's default setting (defaults, by query type). Returns the figure's name and
    the references of the panels (None for the figure as it is), or None if a query type's results or gnuplot script are
    missing or (relative) its default setting was not tested."""
    panels = []
    for query, title in COMBINED_QUERIES:
        base = os.path.join(stats_dir, "%s_parameter_%s_%s_query" % (network, test, query))
        if not os.path.exists(base + ".txt") or not os.path.exists(base + ".gp"):
            return None
        gp_text = open(base + ".gp").read()
        blocks = read_blocks(base + ".txt")
        if len(blocks) != len(re.findall(r" index \d+ title ", gp_text)):
            print("Skipping %s: the blocks of %s.txt do not match the series of its figure" % (test, base))
            return None
        if relative:
            converted = relative_blocks(blocks, test, defaults[query]) if defaults.get(query) else None
            if converted is None:
                print("Skipping the relative %s figure: the %s default setting (par_%s_defaults: %s) is not among the tested values"
                      % (test, query, query, ",".join("%s=%s" % item for item in defaults.get(query, {}).items()) or "not set"))
                return None
            write_blocks(base + "_relative.txt", converted)
            gp_text = gp_text.replace("'%s.txt'" % base, "'%s_relative.txt'" % base)
        lines = [line.strip() for line in gp_text.splitlines()]
        panels.append((title, lines, "%s %s" % (query, reference_label(test, defaults[query])) if relative else None))
    first_lines = panels[0][1]
    style = next(line for line in first_lines if line.startswith("set style line"))
    series = re.findall(r" index \d+ title ('[^']*') with linespoints ls (\d+)", next(l for l in first_lines if l.startswith("plot ")))
    fig_base = "%s_parameter_%s_query_%scombined" % (network, test, "relative_" if relative else "")
    width, height = 5.0 * len(panels), 4.6  # Inches: each panel about as wide as a single figure, plus room for the key
    # Panel edges in screen coordinates: the left margin holds the y label and tics, the bottom the x label and tics, and the top
    # the panel titles and the key, which is drawn by an empty plot across the full width
    left, right, gap, bottom, top, key_bottom = 0.075, 0.975, 0.06, 0.18, 0.76, 0.86
    panel_width = (right - left - gap * (len(panels) - 1)) / len(panels)
    script = ["reset",
              "set terminal postscript enhanced eps 32 color size %g,%g" % (width, height),
              "set output '%s'" % os.path.join(fig_dir, fig_base + ".eps"),
              "set border linewidth 2",
              style,
              "unset key",
              "set multiplot"]
    for i, (title, lines, _) in enumerate(panels):
        panel_left = left + i * (panel_width + gap)
        script.append("set lmargin at screen %.4f; set rmargin at screen %.4f; set bmargin at screen %.4f; set tmargin at screen %.4f"
                      % (panel_left, panel_left + panel_width, bottom, top))
        script.append("set title '%s' font ',34'" % title)
        script.extend(line for line in lines if line.startswith(("set logscale x", "set xlabel", "set xtics", "set autoscale xfix")))
        if relative:
            # The label names no reference, which can be too long for the axis (e.g. Border Minmax); main() prints it for the caption
            script.append("set ylabel 'Relative Query Time'" if i == 0 else "unset ylabel")
            # The default series, flat at 1, is often the lowest: a margin keeps it off the panel's border
            script.append("set offsets 0, 0, graph 0.05, graph 0.05")
        else:
            # Each panel keeps the y axis of its own figure (some are on a log scale), and the first panel's label
            script.extend(line for line in lines if line.startswith(("set logscale y", "set format y")))
            script.append(next(line for line in lines if line.startswith("set ylabel")) if i == 0 else "unset ylabel")
        script.append(re.sub(r"using (\d+):\d+ index", r"using \1:%d index" % TIME_COL,
                             next(line for line in lines if line.startswith("plot "))))
        script.append("unset logscale x; unset logscale y; set format y; set offsets 0, 0, 0, 0")
    # The common key: an empty plot across the top whose entries are the series of the panels
    script += ["set lmargin at screen 0; set rmargin at screen 1; set bmargin at screen %.4f; set tmargin at screen 1" % key_bottom,
               "unset title; unset xlabel; unset ylabel; unset tics; unset border; set xrange [0:1]; set yrange [0:1]",
               "set key center center horizontal maxcols %d samplen 3 spacing 1.1 width 1 font ',28'" % len(series),
               # NaN draws nothing but adds the key entry (a plot of keyentry items only trips an assertion in gnuplot 6.0)
               "plot " + ", ".join("NaN with linespoints ls %s title %s" % (ls, title) for title, ls in series),
               "unset multiplot"]
    gp_out = os.path.join(stats_dir, fig_base + ".gp")
    with open(gp_out, "w") as f:
        f.write("\n".join(script) + "\n")
    subprocess.run(["gnuplot", gp_out], check=True)
    return fig_base, "; ".join(panel[2] for panel in panels) if relative else None


def query_defaults(script_dir, network):
    """Each query type's test defaults in the network's dataset config (par_<query>_defaults), as dictionaries by query type
    (empty for a query type without them, or for every query type if the network has no dataset config)."""
    names = ["par_%s_defaults" % query for query, _ in COMBINED_QUERIES]
    command = '. globalVariables "" --quiet_globalVariables "$1" > /dev/null && ' + " && ".join('echo "${%s:-}"' % name for name in names)
    result = subprocess.run(["bash", "-c", command, "bash", network], cwd=script_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print("No dataset config for %s, so its relative figures are skipped:\n%s" % (network, result.stderr.strip()))
        return {}
    return {query: dict(pair.split("=", 1) for pair in line.split(",") if "=" in pair)
            for (query, _), line in zip(COMBINED_QUERIES, result.stdout.splitlines())}


def main():
    if len(sys.argv) > 1:
        if sys.argv[1] in ("-h", "--help") and len(sys.argv) == 2:
            print(USAGE)
            sys.exit(0)
        print("Unknown argument %s\n\n%s" % (sys.argv[1], USAGE), file=sys.stderr)
        sys.exit(1)
    # The directories of globalVariables, sourced from this script's directory as the other scripts do
    script_dir = os.path.dirname(os.path.abspath(__file__))
    settings = subprocess.run(["bash", "-c", '. globalVariables "" --quiet_globalVariables > /dev/null && echo "$stats_dir" && echo "$fig_dir"'],
                              cwd=script_dir, capture_output=True, text=True, check=True)
    stats_dir, fig_dir = settings.stdout.splitlines()[:2]
    pattern = re.compile(r"^(.+)_parameter_(%s)_(%s)_query\.txt$" % ("|".join(TESTS), "|".join(q for q, _ in COMBINED_QUERIES)))
    networks = sorted({match.group(1) for match in map(pattern.match, os.listdir(stats_dir)) if match})
    for network in networks:
        defaults = query_defaults(script_dir, network)
        for test in TESTS:
            for relative in (False, True):
                result = combined_figure(stats_dir, fig_dir, network, test, relative, defaults)
                if result:
                    print("%s.eps%s" % (result[0], ": relative to the defaults %s" % result[1] if relative else ""))


if __name__ == "__main__":
    main()

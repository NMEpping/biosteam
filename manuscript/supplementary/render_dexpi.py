"""
Draw a DEXPI process model as a block flow diagram with pyDEXPI.

A DEXPI process model describes steps and streams but not necessarily a
drawing. This script lays the steps out from left to right, adds the
drawing to the model as a DEXPI ``Diagram`` (blocks, port positions,
connector lines, flow arrows and labels), and uses pyDEXPI to

* save the model including the diagram as DEXPI XML, and
* render the diagram as SVG (and PDF, if ``cairosvg`` is installed).

Usage::

    python render_dexpi.py <process model>.xml [output stem]

Requires pyDEXPI (https://github.com/process-intelligence-research/pyDEXPI),
which is licensed under the AGPL-3.0. This script is a stand-alone tool
and is not part of BioSTEAM.
"""
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

from pydexpi import load_xml_file
from pydexpi.exports.svg_export import DrawDiagram
from pydexpi.loaders.xml_serializer import DexpiXmlSerializer
from pydexpi.pydexpi_classes import pydexpi_classes as dexpi
from pydexpi.toolkits.process_drawing import (
    BLOCK_SHAPE_PREFIX,
    BLOCK_SYMBOL_REGISTRATION,
)

# Sheet units; y runs down.
BLOCK_WIDTH, BLOCK_HEIGHT = 46., 18.
TERMINAL_WIDTH = 46.
COLUMN_PITCH, ROW_PITCH = 78., 30.
MARGIN = 12.
FONT = 'Latin Modern Roman'
TITLE_SIZE, TEXT_SIZE = 4.2, 3.2

BLACK = dexpi.Color(r=30, g=30, b=30)
WHITE = dexpi.Color(r=255, g=255, b=255)
GRAY = dexpi.Color(r=168, g=150, b=140)        # TUHH warm gray
TURQUOISE = dexpi.Color(r=0, g=193, b=212)     # TUHH turquoise
PURPLE = dexpi.Color(r=114, g=0, b=254)        # TUHH purple
TERMINALS = (dexpi.Source, dexpi.Sink)


def point(x, y):
    return dexpi.Point(x=float(x), y=float(y))


def stroke(color=BLACK, width=0.5):
    return dexpi.Stroke(color=color, dashStyle=dexpi.DashStyle.SOLID,
                        width=width)


def text(content, x, y, size=TEXT_SIZE, color=BLACK,
         alignment=dexpi.TextAlignment.CENTERCENTER):
    return dexpi.Text(
        alignment=alignment, color=color, font=FONT, position=point(x, y),
        rotation=0., size=size, text=content)


def block_shape(width, height, color):
    """Return a block shape: a rectangle centred on its origin."""
    w, h = width / 2., height / 2.
    outline = dexpi.Polygon(
        fillStyle=dexpi.FillStyle.TRANSPARENT,
        points=[point(-w, -h), point(w, -h), point(w, h), point(-w, h)],
        stroke=stroke(color, 0.7))
    return dexpi.Shape(
        id=(f'{BLOCK_SHAPE_PREFIX}_{width:g}x{height:g}'
            f'_{color.r}_{color.g}_{color.b}'),
        name=f'Process block {width:g} x {height:g}',
        symbolRegistrationNumber=BLOCK_SYMBOL_REGISTRATION,
        primitives=[outline])


def arrow_shape():
    """Return an arrowhead pointing in the +x direction, tip at origin."""
    head = dexpi.Polygon(
        fillStyle=dexpi.FillStyle.SOLID,
        points=[point(0., 0.), point(-3., -1.2), point(-3., 1.2)],
        stroke=stroke(BLACK, 0.2))
    return dexpi.Shape(
        id='flow_arrow', name='Flow arrow',
        symbolRegistrationNumber='flow arrow', primitives=[head])


def step_columns(steps, connections):
    """Return the column of each step: the longest path from a feed."""
    owner = {id(port): step for step in steps for port in step.ports}
    downstream = defaultdict(list)
    for connection in connections:
        source = owner.get(id(connection.source))
        target = owner.get(id(connection.target))
        if source is not None and target is not None:
            downstream[id(source)].append(target)
    columns = {}

    def visit(step, column, path):
        if id(step) in path:  # Recycle: keep the column already assigned
            return
        if columns.get(id(step), -1) >= column:
            return
        columns[id(step)] = column
        for target in downstream[id(step)]:
            visit(target, column + 1, path | {id(step)})

    has_inlet = {
        id(step) for step in steps
        if any(port.nominalDirection == dexpi.PortDirection.INLET
               for port in step.ports)}
    for step in steps:
        if id(step) not in has_inlet:
            visit(step, 0, frozenset())
    for step in steps:  # Steps that are only part of a cycle
        columns.setdefault(id(step), 0)
    # Products are drawn in the last column
    last = max(columns.values())
    for step in steps:
        if isinstance(step, dexpi.Sink):
            columns[id(step)] = last
    return columns


def add_diagram(model):
    """Add a block flow diagram of the process model to `model`."""
    process = model.conceptualModel
    steps = list(process.processSteps)
    connections = [i for i in process.processConnections
                   if isinstance(i, dexpi.Stream)]
    columns = step_columns(steps, connections)
    rows = defaultdict(list)
    for step in steps:
        rows[columns[id(step)]].append(step)
    tallest = max(len(i) for i in rows.values())

    shapes = {}

    def shape_for(step):
        terminal = isinstance(step, TERMINALS)
        key = (TERMINAL_WIDTH if terminal else BLOCK_WIDTH,
               GRAY if terminal else PURPLE)
        if key not in shapes:
            shapes[key] = block_shape(key[0], BLOCK_HEIGHT, key[1])
        return shapes[key], key[0]

    arrow = arrow_shape()
    groups, positions = [], {}
    for column, members in sorted(rows.items()):
        offset = (tallest - len(members)) * ROW_PITCH / 2.
        for row, step in enumerate(members):
            x = MARGIN + BLOCK_WIDTH / 2. + column * COLUMN_PITCH
            y = MARGIN + BLOCK_HEIGHT / 2. + offset + row * ROW_PITCH
            shape, width = shape_for(step)
            symbol = dexpi.Static(elements=[dexpi.ShapeUsage(
                shape=shape, isMirrored=False, position=point(x, y),
                rotation=0., scaleX=1., scaleY=1.)])
            lines = [
                text(str(step.identifier).replace('_', ' '), x, y - 3.4,
                     TITLE_SIZE),
                text(type(step).__name__, x, y + 3.2, color=GRAY),
            ]
            nodes = []
            for direction, side in ((dexpi.PortDirection.INLET, -1.),
                                    (dexpi.PortDirection.OUTLET, 1.)):
                ports = [i for i in step.ports
                         if i.nominalDirection == direction]
                for index, port in enumerate(ports):
                    fraction = (index + 1) / (len(ports) + 1)
                    node = dexpi.PortNodePosition(
                        port=port, position=point(
                            x + side * width / 2.,
                            y - BLOCK_HEIGHT / 2. + fraction * BLOCK_HEIGHT))
                    positions[id(port)] = node
                    nodes.append(node)
            groups.append(dexpi.RepresentationGroup(
                represents=step, nodePositions=nodes,
                groups=[symbol, dexpi.Label(elements=lines)]))

    for connection in connections:
        source = positions.get(id(connection.source))
        target = positions.get(id(connection.target))
        if source is None or target is None:
            continue
        x1, y1 = source.position.x, source.position.y
        x2, y2 = target.position.x, target.position.y
        middle = (x1 + x2) / 2.
        inner = []
        if abs(y1 - y2) > 1e-9:
            inner = [point(middle, y1), point(middle, y2)]
        line = dexpi.ConnectorLine(
            source=source, target=target, innerPoints=inner,
            stroke=stroke(TURQUOISE, 0.9))
        head = dexpi.PipeFlowArrow(elements=[dexpi.ShapeUsage(
            shape=arrow, isMirrored=False, position=point(x2, y2),
            rotation=0., scaleX=1., scaleY=1.)])
        label_x = (x1 + middle) / 2. if inner else middle
        label = dexpi.Label(elements=[text(
            str(connection.identifier).replace('_', ' '), label_x, y1 - 1.6,
            alignment=dexpi.TextAlignment.CENTERBOTTOM)])
        groups.append(dexpi.RepresentationGroup(
            represents=connection,
            groups=[dexpi.Static(elements=[line]), head, label]))

    width = 2 * MARGIN + BLOCK_WIDTH + max(rows) * COLUMN_PITCH
    height = 2 * MARGIN + BLOCK_HEIGHT + (tallest - 1) * ROW_PITCH
    model.diagram = dexpi.Diagram(
        represents=process, name='Block flow diagram', groups=groups,
        backgroundColor=WHITE, minX=0., minY=0., maxX=width, maxY=height)
    model.shapeCatalogues = [dexpi.ShapeCatalogue(
        name='Block flow diagram shapes', shapes=[*shapes.values(), arrow])]
    return model.diagram


def render(source, stem=None):
    """Draw the process model in `source`; return the files written."""
    source = Path(source)
    stem = Path(stem) if stem else source.with_name(
        source.name.replace('.dexpi.xml', '').replace('.xml', '') + '_diagram')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = load_xml_file(source, strict=True)
        diagram = add_diagram(model)
        DexpiXmlSerializer().save(model, stem.parent, stem.name + '.dexpi')
        svg = DrawDiagram(diagram, padding=2.).draw_svg(background=False)
    written = [stem.with_name(stem.name + '.dexpi.xml')]
    svg_path = stem.with_suffix('.svg')
    svg_path.write_text(svg, encoding='utf-8')
    written.append(svg_path)
    try:
        import cairosvg
    except ImportError:
        pass
    else:
        pdf_path = stem.with_suffix('.pdf')
        cairosvg.svg2pdf(bytestring=svg.encode('utf-8'),
                         write_to=os.fspath(pdf_path))
        written.append(pdf_path)
    return written


if __name__ == '__main__':
    for path in render(*sys.argv[1:3]):
        print(path)

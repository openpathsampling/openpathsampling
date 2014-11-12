import svgwrite
import sys
import argparse
import os
from storage import Storage

from orderparameter import OP_Function
import mdtraj as md
import numpy as np


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Analyze a file.')
    parser.add_argument('--rejected', dest='rejected', action='store_const',
                   const=True, default=False,
                   help='show also rejected paths')
    parser.add_argument('file', metavar='file.nc', help='an integer for the accumulator')
    parser.add_argument('--show', metavar='orderparameter/snapshot/configuration/momentum', default='',
                   help='show an orderparameter or an id of the snapshot/configuration/momentum')

    parser.add_argument('--pdf', dest='pdf', action='store_const',
                   const=True, default=False,
                   help='create a pdf instead of a svg. Uses wkhtmltopdf to make apdf using webkit and qt')

    parser.add_argument('--state', nargs=2, action='append',
                        help='add a background coloring to s snapshot using the given orderparameter.',
                        metavar=('orderparameter', 'color'))

    parser.add_argument('--show-degree', dest='in_degree', action='store_const',
                   const=True, default=False,
                   help='shows orderparameters in degree and not plain')

    args = parser.parse_args()

    rejected = args.rejected
    file = args.file

    degrees = 180/3.14159 # psi reports in radians; I think in degrees

    if not os.path.isfile(file):
        print file, 'does not exist ! ENDING!'
        exit()

    storage = Storage(
         filename = file,
         mode = 'a'
         )


    if args.show == 'snapshot':
        def show_id(s):
            return s.idx[storage]

        show_op = show_id
    elif args.show == 'configuration':
        def show_id(s):
            return s.configuration.idx[storage]

        show_op = show_id
    elif args.show == 'momentum':
        def show_id(s):
            return s.momentum.idx[storage]

        show_op = show_id
    elif args.show == '':
        def show_id(s):
            return ''

        show_op = show_id
    else:
        show_op = storage.cv.load(args.show)

    if args.in_degree:
        inner = show_op
        def func(s):
            return degrees * inner(s)

        show_op = func

    matrix = [[]]

    def block(x,y, color = "blue", text = "", snapshot = None):
        global matrix
        svg_document.add(svg_document.rect(
            insert = (start_x + (x +0.05) * scale_x ,start_y + (y - 0.3) * scale_y),
            size = (0.9 * scale_x, 0.6 * scale_y),
            fill = color,
        ))
        svg_document.add(svg_document.circle(
                                              center = (start_x + x * scale_x,start_y + y * scale_y),
                                                  r = 0.05 * scale_x,
                                                 stroke_width = "0",
                                                  stroke = color,
                                                  fill = color
                                                  ))
        svg_document.add(svg_document.circle(
                                              center = (start_x + (x + 1) * scale_x,start_y + y * scale_y),
                                                  r = 0.05 * scale_x,
                                                 stroke_width = "0",
                                                  stroke = color,
                                                  fill = color
                                                  ))
        svg_document.add(svg_document.text(
                                        text = ("     " + (str(text)[:5]))[-5:],
                                        insert = (start_x + (x + 0.5) * scale_x,start_y + (y + 0.05) * scale_y),
                                        text_anchor = 'middle',
                                        font_size = 0.3*scale_y,
                                        alignment_baseline = 'middle',
                                        font_family = 'Futura',
                                        fill = 'white'
                                        ))

    p_x = dict()
    p_y = dict()

    # First run to determine plot size

    t_count = 0

    min_x = 0
    max_x = 0

    matrix = [[None] * 1000 for n in range(500)]

    for o_idx in range(0, storage.sample.count()):
        sample = storage.sample.load(o_idx)
        if hasattr(sample.details, 'start_point'):
            length = len(sample.details.final)

            old_traj = sample.details.start_point.trajectory
            old_index = sample.details.start_point.index
            old_conf = old_traj[old_index].configuration

            new_traj = sample.details.final_point.trajectory
            new_index = sample.details.final_point.index
            new_conf = new_traj[new_index].configuration

            if sample.trajectory is new_traj or rejected:
                accepted = sample.trajectory is new_traj
                t_count += 0.8
                if not old_conf in p_x:
                    for pos, snapshot in enumerate(old_traj):
                        conf = snapshot.configuration
                        p_x[conf] = pos
                        p_y[conf] = t_count

                        pos_x = p_x[conf]
                        pos_y = p_y[conf]

                        x = int(pos_x)
                        y = int(pos_y / 0.8 + 0.4)
                        matrix[y][x + 500] = snapshot


                    t_count += 0.8

                shift = p_x[old_conf] - new_index

                min_x = min(min_x, shift)
                max_x = max(max_x, shift + len(new_traj) - 1)

                for pos, snapshot in enumerate(new_traj):
                    conf = snapshot.configuration
                    if not conf in p_y:
                        p_y[conf] = t_count
                        p_x[conf] = shift + pos

                        pos_x = p_x[conf]
                        pos_y = p_y[conf]

                        x = int(pos_x)
                        y = int(pos_y / 0.8 + 0.4)
                        matrix[y][x + 500] = snapshot


        else:
            accepted = sample.trajectory is new_traj
            t_count += 2 * 0.8
            # Show full trajectory once
            for pos, snapshot in enumerate(sample.trajectory):
                conf = snapshot.configuration
                if conf not in p_x:
                    p_x[conf] = pos
                p_y[conf] = t_count

                pos_x = p_x[conf]
                pos_y = p_y[conf]

                x = int(pos_x)
                y = int(pos_y / 0.8 + 0.4)
                matrix[y][x + 500] = snapshot


            shift = p_x[old_conf] - new_index

    y_max = int(t_count / 0.8) + 1


    svg_document = svgwrite.Drawing(
        filename = "tree.svg"
    )

    p_x = dict()
    p_y = dict()

    scale_x = 18
    scale_y = 18

    totalwidth = (max_x - min_x + 5) * (scale_x)
    left = (-min_x + 2.5) * (scale_x)

    svg_document['width'] = str(totalwidth)+'px'

    start_x = left
    start_y = 0

    t_count = 0

    lightcolor = "lightgray"

    op_names = {}

    if args.state is not None:
        op_names = { arg[0] : arg[1] for arg in args.state }

    ops = {op : storage.cv.load(op) for op in op_names.keys() }

    for y in range(0, y_max):
        rr = { op_name : None for op_name in op_names.keys() }
        for x in range(0, (max_x - min_x + 1)):
            yp = y * 0.8
            xp = x + min_x
            for r in rr:
                op = ops[r]
                if matrix[y][xp + 500] is not None and bool(op(matrix[y][xp + 500])):
                    if rr[r] is None:
                        rr[r] = xp
                else:
                    if rr[r] is not None:
                        svg_document.add(svg_document.rect(
                            insert = (start_x + rr[r] * scale_x ,start_y + (yp - 0.35) * scale_y),
                            size = ((1.0 + xp - rr[r] - 1) * scale_x, 0.7 * scale_y),
                            fill = op_names[r],
                            stroke = op_names[r],
#                            stroke_dasharray='1,0.2',
                            stroke_width = 0.8
                        ))
                        rr[r] = None

        for r in rr:
            if rr[r] is not None:
                svg_document.add(svg_document.rect(
                    insert = (start_x + (rr[r]) * scale_x ,start_y + (yp - 0.35) * scale_y),
                    size = ((1.0 + xp - rr[r]) * scale_x, 0.7 * scale_y),
                    fill = 'none',
                    stroke = op_names[r],
#                    stroke_dasharray='1.0,0.2',
                    stroke_width = 0.8
                ))


    for o_idx in range(0, storage.sample.count()):
        sample = storage.sample.load(o_idx)
        if hasattr(sample.details, 'start_point'):
            length = len(sample.details.final)

            old_traj = sample.details.start_point.trajectory
            old_index = sample.details.start_point.index
            old_conf = old_traj[old_index].configuration

            new_traj = sample.details.final_point.trajectory
            new_index = sample.details.final_point.index
            new_conf = new_traj[new_index].configuration

            if sample.trajectory is new_traj or rejected:
                accepted = sample.trajectory is new_traj
                t_count += 0.8
                if not old_conf in p_x:
                    for pos, snapshot in enumerate(old_traj):
                        conf = snapshot.configuration
                        p_x[conf] = pos
                        p_y[conf] = t_count

                        pos_x = p_x[conf]
                        pos_y = p_y[conf]

                        block(pos_x, pos_y, "black", show_op(snapshot))

                    t_count += 0.8

                shift = p_x[old_conf] - new_index

                fontcolor = "black"

                if sample.details.mover.name == "BackwardShootMover":
                    color = "green"
                    if not accepted:
                        color = lightcolor
                        fontcolor = lightcolor
                    svg_document.add(svg_document.line(
                        start = (start_x + (shift + new_index) * scale_x + 0,start_y + (p_y[old_conf] + 0.05) * scale_y),
                        end = (start_x  + (shift + new_index) * scale_x + 0,start_y + (t_count - 0.05) * scale_y),
                        stroke_width = 0.05 * scale_x,
                        stroke = color,
                    ))
                    svg_document.add(svg_document.text(
                                            text = str(new_traj.idx[storage]) + 'b',
                                            insert = (start_x + (shift + 0 - 0.2) * scale_x,start_y + (t_count + 0.05) * scale_y),
                                            text_anchor = 'end',
                                            font_size = 0.4*scale_y,
                                            alignment_baseline = 'middle',
                                            font_family = 'Futura',
                                            fill = fontcolor
                                            ))
                else:
                    color = "red"
                    if not accepted:
                        color = lightcolor
                        fontcolor = lightcolor

                    svg_document.add(svg_document.line(
                        start = (start_x + (shift + new_index + 1) * scale_x,start_y + 0.5 + p_y[old_conf] * scale_y),
                        end = (start_x  + (shift + new_index + 1) * scale_x,start_y - 0.5 + t_count * scale_y),
                        stroke_width = 0.05 * scale_x,
                        stroke = color,
                    ))
                    svg_document.add(svg_document.text(
                                            text = str(new_traj.idx[storage]) + 'f',
                                            insert = (start_x + (shift + len(new_traj) + 0.2) * scale_x,start_y + (t_count + 0.05) * scale_y),
                                            text_anchor = 'start',
                                            font_size = 0.4*scale_y,
                                            alignment_baseline = 'middle',
                                            font_family = 'Futura',
                                            fill = fontcolor
                                            ))

                if not  accepted:
                    color = lightcolor

                for pos, snapshot in enumerate(new_traj):
                    conf = snapshot.configuration
                    if not conf in p_y:
                        p_y[conf] = t_count
                        p_x[conf] = shift + pos

                        pos_x = p_x[conf]
                        pos_y = p_y[conf]

                        block(pos_x, pos_y, color, show_op(snapshot))

        else:
            accepted = sample.trajectory is new_traj
            t_count += 2 * 0.8
            # Show full trajectory once
            for pos, snapshot in enumerate(sample.trajectory):
                conf = snapshot.configuration
                if conf not in p_x:
                    p_x[conf] = pos
                p_y[conf] = t_count

                pos_x = p_x[conf]
                pos_y = p_y[conf]

                block(pos_x, pos_y, color, show_op(snapshot))

            shift = p_x[old_conf] - new_index

    svg_document['height'] = str(int(start_y + (t_count + 2)*scale_y))+'px'

    svg_document.save()

    html_file = '<!DOCTYPE html><html style="margin:0px; padding:0px;">' + svg_document.tostring() + '<body style="margin:0px; padding:0px;"></body></html>'

    h = int(start_y + (t_count + 2)*scale_y)
    w = int(totalwidth)

    h_margin = 20.0
    v_margin = 22.0

    page_height =str((210.0 - h_margin) * h/w + v_margin)+'mm'

    with open('tree.html', 'w') as f:
        f.write(html_file)
    f.closed

    bashCommand = "open -a Safari tree.svg"
#    os.system(bashCommand)

    if args.pdf:
        bashCommand = "open -a Safari tree.html"
    #    os.system(bashCommand)

        bashCommand = "wkhtmltopdf -l --page-height " + page_height + " --page-width 210.0mm tree.html tree.pdf"
        os.system(bashCommand)

        bashCommand = "open tree.pdf"
        os.system(bashCommand)

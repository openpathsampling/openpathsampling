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
    parser.add_argument('--psi', dest='psi', action='store_const',
                   const=True, default=False,
                   help='show psi angle instead of ID')
#    parser.add_argument('--phi', dest='phi', action='store_const',
#                   const=True, default=False,
#                   help='show phi angle instead of ID')

    parser.add_argument('--pdf', dest='rejected', action='store_const',
                   const=True, default=False,
                   help='create a pdf instead of a svg. Uses wkhtmltopdf to make apdf using webkit and qt')


    parser.add_argument('--state', nargs=2, action='append')

    args = parser.parse_args()

    rejected = args.rejected
    file = args.file
    show_psi = args.psi

    if not os.path.isfile(file):
        print file, 'does not exist ! ENDING!'
        exit()

    storage = Storage(
                                                     filename = file,
                                                     mode = 'a'
                                                     )

    psi_atoms = [6,8,14,16]
    psi = OP_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
      indices=[psi_atoms])


    def headline(s):
        print
        print "###################################################################################"
        print "##", s.upper()
        print "###################################################################################"
        print

    def line(a, b):
        print '    {:<32} : {:<30}'.format(a,b)

    def nline(n, a, b):
        print '     {:>4}] {:<25} : {:<30}'.format(n,a,b)

    def print_traj(name, traj_obj):
        traj = storage.trajectory.configuration_indices(traj_obj.idx[storage])
        sys.stdout.write("      {:>10}:  {:>5} frames [".format(name, str(len(traj))))
        old_idx = -2
        count = 0
        for idx in traj:
            if idx == old_idx + 1 or idx == old_idx - 1:
                count += 1
            else:
                if count > 1:
                    sys.stdout.write(" <" + str(count - 1) + ">")
                if old_idx >= 0 and count > 0:
                    sys.stdout.write(" " + str(old_idx))
                sys.stdout.write(" " + str(idx))
                count = 0
            old_idx = idx

        if count > 1:
            sys.stdout.write(" ... <" + str(count - 1) + "> ...")
        if count > 0:
            sys.stdout.write(" " + str(old_idx))


        sys.stdout.write(" ]\n")

    headline("General")

    line("Filename", file)
    line("Size", str(os.path.getsize(file) / 1024 / 1024) + " MB")

    matrix = [[]]

    def block(x,y, color = "blue", idx = "", snapshot = None):
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
                                        text = str(idx),
                                        insert = (start_x + (x + 0.5) * scale_x,start_y + (y + 0.05) * scale_y),
                                        text_anchor = 'middle',
                                        font_size = 0.4*scale_y,
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

    degrees = 180/3.14159 # psi reports in radians; I think in degrees

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

                        if show_psi:
                            block(pos_x, pos_y, "black", str(int((degrees * psi(conf))) % 360 ), snapshot=snapshot)
                        else:
                            block(pos_x, pos_y, "black", conf.idx[storage], snapshot=snapshot)

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

                        if show_psi:
                            block(pos_x, pos_y, color, str(int((degrees * psi(conf))) % 360 ), snapshot=snapshot)
                        else:
                            block(pos_x, pos_y, color, conf.idx[storage], snapshot=snapshot)

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

                if show_psi:
                    block(pos_x, pos_y, "black", str(int((degrees * psi(conf))) % 360 ), snapshot=snapshot)
                else:
                    block(pos_x, pos_y, "black", conf.idx[storage], snapshot=snapshot)

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

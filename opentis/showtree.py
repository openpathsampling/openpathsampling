import svgwrite
import sys
import argparse
import os
from storage import Storage

if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Analyze a file.')
    parser.add_argument('file', metavar='file.nc', help='an integer for the accumulator')

    args = parser.parse_args()
    file = args.file

    if not os.path.isfile(file):
        print file, 'does not exist ! ENDING!'
        exit()

    storage = Storage(
                                                     filename = file,
                                                     mode = 'a'
                                                     )

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

    def block(x,y, color = "blue", idx = ""):
        svg_document.add(svg_document.line(
            start = (start_x + x * scale + 0.5,start_y + y * scale),
            end = (start_x  + (x + 1) * scale - 0.5,start_y + y * scale),
            stroke_width = "5",
            stroke = color,
        ))
        svg_document.add(svg_document.circle(
                                              center = (start_x + x * scale,start_y + y * scale),
                                                  r = 0.5,
                                                 stroke_width = "0",
                                                  stroke = color,
                                                  fill = color
                                                  ))
        svg_document.add(svg_document.circle(
                                              center = (start_x + (x + 1) * scale,start_y + y * scale),
                                                  r = 0.5,
                                                 stroke_width = "0",
                                                  stroke = color,
                                                  fill = color
                                                  ))
        svg_document.add(svg_document.text(
                                        text = str(idx),
                                        insert = (start_x + (x + 0.5) * scale,start_y + (y + 0.05) * scale),
                                        text_anchor = 'middle',
                                        font_size = '4',
                                        alignment_baseline = 'middle',
                                        font_family = 'Futura',
                                        fill = 'white'
                                        ))

    svg_document = svgwrite.Drawing(
        filename = "test-svgwrite.svg"
    )

    p_x = dict()
    p_y = dict()

    start_x = 500
    start_y = 100

    scale = 10
    t_count = 0

    lightcolor = "lightgray"

    for o_idx in range(1, storage.sample.count() + 1):
        sample = storage.sample.load(o_idx)
        length = len(sample.details.final)

        old_traj = sample.details.start_point.trajectory
        old_index = sample.details.start_point.index
        old_conf = old_traj[old_index].configuration

        new_traj = sample.details.final_point.trajectory
        new_index = sample.details.final_point.index
        new_conf = new_traj[new_index].configuration

        if sample.trajectory is new_traj or True:
            accepted = sample.trajectory is new_traj
            t_count += 1
            if not old_conf in p_x:
                for pos, snapshot in enumerate(old_traj):
                    conf = snapshot.configuration
                    p_x[conf] = pos
                    p_y[conf] = t_count

                    pos_x = p_x[conf]
                    pos_y = p_y[conf]

                    block(pos_x, pos_y, "blue", conf.idx[storage])

                t_count += 1

            shift = p_x[old_conf] - new_index

            if sample.details.mover.name == "BackwardShootMover":
                color = "green"
                if not accepted:
                    color = lightcolor
                svg_document.add(svg_document.line(
                    start = (start_x + (shift + new_index) * scale + 0,start_y + 0.5 + p_y[old_conf] * scale),
                    end = (start_x  + (shift + new_index) * scale + 0,start_y - 0.5 + t_count * scale),
                    stroke_width = "0.5",
                    stroke = color,
                ))
                svg_document.add(svg_document.text(
                                        text = str(new_traj.idx[storage]) + 'b',
                                        insert = (start_x + (shift + 0 - 0.2) * scale,start_y + (t_count + 0.05) * scale),
                                        text_anchor = 'end',
                                        font_size = '4',
                                        alignment_baseline = 'middle',
                                        font_family = 'Futura',
                                        fill = 'black'
                                        ))

            else:
                color = "red"
                if not accepted:
                    color = lightcolor

                svg_document.add(svg_document.line(
                    start = (start_x + (shift + new_index + 1) * scale,start_y + 0.5 + p_y[old_conf] * scale),
                    end = (start_x  + (shift + new_index + 1) * scale,start_y - 0.5 + t_count * scale),
                    stroke_width = "0.5",
                    stroke = color,
                ))
                svg_document.add(svg_document.text(
                                        text = str(new_traj.idx[storage]) + 'f',
                                        insert = (start_x + (shift + len(new_traj) + 0.2) * scale,start_y + (t_count + 0.05) * scale),
                                        text_anchor = 'start',
                                        font_size = '4',
                                        alignment_baseline = 'middle',
                                        font_family = 'Futura',
                                        fill = 'black'
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

                    block(pos_x, pos_y, color, conf.idx[storage])

    svg_document.save()
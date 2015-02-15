import svgwrite
import sys
import argparse
import os
from openpathsampling.storage import Storage

from openpathsampling.orderparameter import OP_Function
import mdtraj as md
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

    headline("General")

    line("Filename", file)
    line("Size", str(os.path.getsize(file) / 1024 / 1024) + " MB")

    def block(x,y, color = "blue", idx = ""):
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

    svg_document = svgwrite.Drawing(
        filename = "tree.svg"
    )

    p_x = dict()
    p_y = dict()

    scale_x = 18
    scale_y = 18

    totalwidth = 1800
    left = 0

    svg_document['width'] = str(totalwidth)+'px'

    start_x = left
    start_y = 50

    t_count = 0

    lightcolor = "lightgray"

    degrees = 180/3.14159 # psi reports in radians; I think in degrees

    ensembles = map(storage.ensemble.load, range(storage.ensemble.count()))
    ensembles_y = {ensemble : idx for idx, ensemble in enumerate(ensembles) if ensemble.name != ''}

    left = 100

    state = dict()
    changed = set()

    for ensemble, pos_y in ensembles_y.iteritems():
        svg_document.add(svg_document.text(
                text = '[' + str(ensemble.idx[storage]) + '] ' + str(ensemble.name) ,
                insert = (start_x + (0 + 0.5) * scale_x, start_y + (pos_y + 0.05) * scale_y),
                text_anchor = 'start',
                font_size = 0.6*scale_y,
                alignment_baseline = 'middle',
                font_family = 'Futura',
                fill = 'black'
                ))

        state[ensemble] = None

    start_x = left

    svg_document.add(svg_document.text(
            text = 'time' ,
            insert = (start_x + (0 + 0.5) * scale_x, start_y + (-1 + 0.05) * scale_y),
            text_anchor = 'end',
            font_size = 0.6*scale_y,
            alignment_baseline = 'middle',
            font_family = 'Futura',
            fill = 'black'
            ))


    last_x = -1


    for o_idx in range(0, storage.sample.count()):

        sample = storage.sample.load(o_idx)
        ensemble = sample.ensemble
        length = len(sample.details.final)

        p_y[sample] = ensembles_y[ensemble]
        p_x[sample] = sample.step

        pos_x = p_x[sample]
        pos_y = p_y[sample]

        if sample.details.mover.name == 'BootstrapEnsembleChangeMove':
            svg_document.add(svg_document.line(
                    start = (start_x + (pos_x) * scale_x + 0,start_y + (pos_y - 1 + 0.05) * scale_y),
                    end = (start_x  + (pos_x) * scale_x + 0,start_y + (pos_y - 0.05) * scale_y),
                    stroke_width = 0.05 * scale_x,
                    stroke = "black",
                ))

            svg_document.add(svg_document.line(
                    start = (start_x + (pos_x - 0.1) * scale_x + 0,start_y + (pos_y - 0.05 - 0.5) * scale_y),
                    end = (start_x  + (pos_x) * scale_x + 0,start_y + (pos_y + 0.05 - 0.5) * scale_y),
                    stroke_width = 0.05 * scale_x,
                    stroke = "black",
                ))
            svg_document.add(svg_document.line(
                    start = (start_x + (pos_x + 0.1) * scale_x + 0,start_y + (pos_y - 0.05 - 0.5) * scale_y),
                    end = (start_x  + (pos_x) * scale_x + 0,start_y + (pos_y + 0.05 - 0.5) * scale_y),
                    stroke_width = 0.05 * scale_x,
                    stroke = "black",
                ))

        if sample.details.mover.name == 'ReplicaExchange':
            middle_y = 0.5 * (ensembles_y[sample.details.ensembles[0]] + ensembles_y[sample.details.ensembles[1]])
            svg_document.add(svg_document.line(
                    start = (start_x + (pos_x) * scale_x + 0,start_y + (middle_y - 0.5 + 0.05) * scale_y),
                    end = (start_x  + (pos_x) * scale_x + 0,start_y + (middle_y + 0.5 - 0.05) * scale_y),
                    stroke_width = 0.05 * scale_x,
                    stroke = "black",
                ))
            svg_document.add(svg_document.circle(
                  center = (start_x + (pos_x) * scale_x,start_y + (middle_y) * scale_y),
                      r = 0.1 * scale_x,
                     stroke_width = "0",
                      stroke = "black",
                      fill = "black"
                      ))

        if pos_x > last_x:
            for ens, pos_yy in ensembles_y.iteritems():
                svg_document.add(svg_document.line(
                    start = (start_x + (pos_x) * scale_x + 0,start_y + (pos_yy) * scale_y),
                    end = (start_x  + (pos_x + 1) * scale_x + 0,start_y + (pos_yy) * scale_y),
                    stroke_width = 0.025 * scale_x,
                    stroke = "black",
                    stroke_dasharray=0.1*scale_x
                ))

            if last_x >=0:
                for ens, pos_yy in ensembles_y.iteritems():
                    if ens not in changed and state[ens] is not None:
                        block(last_x, pos_yy, lightcolor, str(int(state[ens].trajectory.idx[storage])))

            changed = set()

            last_x = pos_x

            block(pos_x, -1, "red", str(int(sample.step)))

        block(pos_x, pos_y, "black", str(int(sample.trajectory.idx[storage])))
        changed.add(ensemble)
        state[ensemble] = sample

    if last_x >=0:
        for ens, pos_yy in ensembles_y.iteritems():
            if ens not in changed and state[ens] is not None:
                block(last_x, pos_yy, lightcolor, str(int(state[ens].trajectory.idx[storage])))


    svg_document.save()

    bashCommand = "open -a Safari tree.svg"
    os.system(bashCommand)

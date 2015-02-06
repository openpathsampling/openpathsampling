import argparse
import os
from openpathsampling.storage import Storage
from openpathsampling.visualize import PathTreeBuilder

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

    parser.add_argument('--html', dest='html', action='store_const',
                   const=True, default=False,
                   help='create a svg embedded in a html')

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
        show_op = storage.snapshot.op_idx
    elif args.show == 'configuration':
        show_op = storage.snapshot.op_configuration_idx
    elif args.show == 'momentum':
        show_op = storage.snapshot.op_momentum_idx
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

    tree = PathTreeBuilder(storage, op=show_op, states=args.state)
    samples = storage.sample.by_ensemble(storage.ensemble.load(4))
    tree.from_samples(samples)

    if args.pdf:
        tree.renderer.save_pdf()

        bashCommand = "open tree.pdf"
        os.system(bashCommand)
    elif args.html:
        tree.renderer.write_html()

        bashCommand = "open -a Safari tree.html"
        os.system(bashCommand)
    else:
        bashCommand = "open -a Safari tree.svg"
    #    os.system(bashCommand)

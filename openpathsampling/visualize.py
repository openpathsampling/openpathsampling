import svgwrite as svg
from svgwrite.container import Group
import os
import openpathsampling as paths

import json
from collections import namedtuple


class TreeRenderer(svg.Drawing):
    def __init__(self, css_style=''):
        super(TreeRenderer, self).__init__()
        self.scale_x = 20.0
        self.scale_y = 20.0

        # Add the CSS Stylesheet
        self.css_style = css_style
        self.defs.add(self.style(
            self.css_style
        ))
        self.horizontal_gap = 0.05

    @staticmethod
    def c(cls):
        return ' '.join(cls)

    def x(self, x):
        return self.w(x)

    def y(self, y):
        return self.h(y)

    def w(self, y):
        return self.scale_x * y

    def h(self, y):
        return self.scale_y * y

    def xy(self, x, y):
        return self.x(x), self.y(y)

    def wh(self, w, h):
        return self.w(w), self.h(h)

    def connector(self, x, y, text="", cls=None):
        if cls is None:
            cls = list()

        cls += ['connector']

        return self.block(x, y, text, False, False, True, True, cls=cls)

    def block(self, x, y, text="",
              extend_right=True, extend_left=True,
              extend_top=False, extend_bottom=False,
              w=1.0, cls=None, data=None):

        if cls is None:
            cls = list()

        cls += ['block']

        padding = self.horizontal_gap

        group = self.g(
            class_=self.c(cls)
        )

        if data is not None:
            group.set_desc(desc=json.dumps(data))

        group.add(self.rect(
            insert=self.xy(x - 0.5 + padding, y - 0.3),
            size=self.wh(1.0 * w - 2 * padding, 0.6),
        ))

        if extend_left:
            group.add(self.circle(
                center=self.xy(x - 0.5, y),
                r=self.w(padding)
            ))
        if extend_right:
            group.add(self.circle(
                center=(self.xy(x + w - 0.5, y)),
                r=self.w(padding)
            ))

        if extend_top:
            group.add(self.circle(
                center=self.xy(x, y - 0.3),
                r=self.w(padding)
            ))
        if extend_bottom:
            group.add(self.circle(
                center=(self.xy(x + w - 1.0, y + 0.3)),
                r=self.w(padding)
            ))

        group.add(self.text(
            text=str(text)[:4],
            insert=self.xy(x + (w - 1.0) / 2.0, y)
        ))

        return group

    def horizontal_region(self, x, y, w=1.0, text="",
                          extend_right=True, extend_left=True, cls=None):

        if cls is None:
            cls = list()

        cls += ['h-region']

        if w == 0:
            return []

        padding = self.horizontal_gap

        group = Group(
            class_=self.c(cls)
        )

        group.add(self.line(
            start=self.xy(x - 0.5 + padding, y),
            end=self.xy(x - 0.5 + w - padding, y)
        ))

        if extend_left:
            group.add(self.circle(
                center=self.xy(x - 0.5, y),
                r=self.w(padding)
            ))
            group.add(self.line(
                start=self.xy(x - 0.5, y - 0.3),
                end=self.xy(x - 0.5, y + 0.3)
            ))

        if extend_right:
            group.add(self.circle(
                center=(self.xy(x + w - 0.5, y)),
                r=self.w(padding)
            ))
            group.add(self.line(
                start=self.xy(x + w - 0.5, y - 0.3),
                end=self.xy(x + w - 0.5, y + 0.3)
            ))

        group.add(self.text(
            text=str(text),
            insert=self.xy(x + (w - 1.0) / 2.0, y),
            class_='shadow'
        ))
        group.add(self.text(
            text=str(text),
            insert=self.xy(x + (w - 1.0) / 2.0, y)
        ))

        return group

    def vertical_region(self, x, y, w=1.0, text="", extend_top=True, extend_bottom=True, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-region']

        padding = self.horizontal_gap
        width = 0.2
        gap = 0.0

        group = Group(
            class_=self.c(cls)
        )

        group.add(self.line(
            start=self.xy(x, y - 0.5 + gap),
            end=self.xy(x, y + w - 1 + 0.5 - gap)
        ))

        if extend_top:
            group.add(self.circle(
                center=self.xy(x, y - 0.5 + gap),
                r=self.w(padding)
            ))
            group.add(self.line(
                start=self.xy(x - 0 * width, y - 0.5 + gap),
                end=self.xy(x + width, y - 0.5 + gap)
            ))

        if extend_bottom:
            group.add(self.circle(
                center=(self.xy(x, y + (w - 1.0) + 0.5 - gap)),
                r=self.w(padding)
            ))
            group.add(self.line(
                start=self.xy(x - 0 * width, y + w - 1.0 + 0.5 - gap),
                end=self.xy(x + width, y + w - 1.0 + 0.5 - gap)
            ))

        group.add(self.text(
            text=str(text),
            insert=self.xy(x - width, y + (w - 1.0) / 2.0)
        ))

        return group

    def shade(self, x, y, w, cls=None, color=None):
        if cls is None:
            cls = list()

        cls += ['shade']

        if color is None:
            return self.rect(
                class_=self.c(cls),
                insert=self.xy(x - 0.5, y + 0.35),
                size=self.wh(w, 0.1)
            )
        else:
            return self.rect(
                class_=self.c(cls),
                insert=self.xy(x - 0.5, y + 0.35),
                size=self.wh(w, 0.1),
                fill=color
            )

    def vertical_connector(self, x, y1, y2, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-connector']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self.xy(x - 0.5, y1 + padding),
            end=self.xy(x - 0.5, y2 - padding)
        )

    def vertical_hook(self, x1, y1, x2, y2, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-hook']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self.xy(x1, y1 + padding + 0.3),
            end=self.xy(x2, y2 - padding - 0.3)
        )

    def horizontal_connector(self, x1, x2, y, cls=None):
        if cls is None:
            cls = list()

        cls += ['h-connector']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self.xy(x1 + 0.5 + padding, y),
            end=self.xy(x2 - 0.5, y)
        )

    def label(self, x, y, text, cls=None):
        if cls is None:
            cls = list()

        cls += ['label']

        group = self.g(
            class_=self.c(cls)
        )

        group.translate(self.x(x), self.y(y))

        group2 = self.g(
            class_='shift'
        )

        group2.add(
            self.text(
                text=str(text),
                insert=(0, 0)
            )
        )

        group.add(
            group2
        )

        return group

    def vertical_label(self, x, y, text, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-label']

        group = self.g(
            class_=self.c(cls)
        )

        group.translate(x, y)

        group.add(
            self.text(
                text=str(text),
                insert=(0, 0),
            )
        )

        return group

    def rectangle(self, x, y, w, h, cls=None):
        if cls is None:
            cls = list()

        return self.rect(
            class_=self.c(cls),
            insert=self.xy(x, y),
            size=self.wh(w, h),
        )

    def to_svg(self):
        return self.tostring()

    def to_html(self):
        svg_source = self.to_svg()
        html = '<!DOCTYPE html><html style="margin:0px; padding:0px; width:100%;">' + \
               svg_source + '<body style="margin:0px; padding:0px;"></body></html>'

        return html

    def _height(self):
        return self.h(self.height) + self.margin * 2

    def _width(self):
        return self.w(self.width) + self.margin * 2

    def write_html(self, file_name='tree.html'):
        with open(file_name, 'w') as f:
            f.write(self.to_html())

    def save_pdf(self, file_name='tree.pdf'):
        with open('tree_xxx.html', 'w') as f:
            f.write(self.to_html())

        bash_command = "open tree_xxx.html " + file_name

        os.system(bash_command)


class Builder(object):
    def __init__(self, additional_option_categories=None):
        options = ['analysis', 'css', 'ui', 'format']
        if additional_option_categories is not None:
            options += additional_option_categories

        option_tuple_class = namedtuple(
            'optionstuple',
            ' '.join(options)
        )

        self.options = option_tuple_class(**{opt: {} for opt in options})

    def svg(self):
        return self.render().tostring()

    def html(self):
        return self.render().tostring()

    def render(self):
        raise NotImplemented('This is a stub class. Use a derived instance!')


class MoveTreeBuilder(Builder):
    def __init__(self):
        super(MoveTreeBuilder, self).__init__()

        self.rejected = False
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()

        self.ensembles = []
        self.pathmover = None

        self.traj_ens_x = dict()
        self.traj_ens_y = dict()

        self.traj_repl_x = dict()
        self.traj_repl_y = dict()

        self.ens_x = list()
        self.repl_x = list()

        self.css_style = vis_css
        self.options.analysis['only_canonical'] = True

        self.doc = None

    def set_ensembles(self, ensembles):
        self.ensembles = ensembles

    def set_mover(self, pathmover):
        self.pathmover = pathmover

    def render(self):
        doc = TreeRenderer(self.css_style)
        self.doc = doc

        level_y = dict()

        self.ens_x = [None] * len(self.ensembles)
        self.repl_x = [None] * len(self.ensembles)

        path = self.pathmover

        group = doc.g(
            class_='tree'
        )

        tree = path.depth_pre_order(lambda this: this, only_canonical=self.options.analysis['only_canonical'])
        total = len(tree)

        for yp, (level, sub_mp) in enumerate(tree):
            x_pos = - level

            sub_type = sub_mp.__class__
            sub_name = sub_type.__name__[:-5]

            if sub_type is paths.SamplePathMoveChange:
                group.add(
                    doc.block(level, yp))

                group.add(
                    doc.label(
                        x_pos,
                        yp,
                        sub_name,
                        cls=['name'] + [sub_type.__name__]
                    )
                )
            else:
                group.add(
                    doc.block(
                        x_pos,
                        yp,
                    )
                )
                group.add(
                    doc.label(
                        x_pos,
                        yp,
                        sub_name
                    )
                )

            if level - 1 in level_y \
                    and level_y[level - 1] == yp - 1:
                group.add(
                    doc.vertical_connector(
                        x_pos + 1,
                        yp,
                        yp - 1
                    )
                )

            if level + 1 in level_y:
                del level_y[level + 1]

            if level in level_y and level_y[level]:
                group.add(
                    doc.vertical_connector(
                        x_pos + 1,
                        yp,
                        level_y[level]
                    )
                )

            level_y[level] = yp

        doc.add(group)

        group = doc.g(
            class_='ensembles'
        )

        for ens_idx, ens in enumerate(self.ensembles):
            txt = chr(ens_idx + 65)

            label = ens.name if hasattr(ens, 'name') else ens.__class__.__name__[:-8]

            group.add(
                doc.label(
                    ens_idx,
                    -1,
                    '[' + txt + '] ' + label,
                    cls=['head']
                )
            )
            group.add(
                doc.vertical_hook(
                    ens_idx,
                    -1,
                    ens_idx,
                    total
                )
            )

        max_level = 0

        for yp, (level, sub_mp) in enumerate(
                path.depth_pre_order(lambda this: this, only_canonical=self.options.analysis['only_canonical'])):
            if level > max_level:
                max_level = level

            in_ens = sub_mp.input_ensembles
            out_ens = sub_mp.output_ensembles

            for ens_idx, ens in enumerate(self.ensembles):
                txt = chr(ens_idx + 65)
                show = False

                if in_ens is None or None in in_ens or ens in in_ens:
                    group.add(
                        doc.connector(
                            ens_idx,
                            yp - 0.15,
                            cls=['input']
                        )
                    )
                    show = True

                if out_ens is None or None in out_ens or ens in out_ens:
                    group.add(
                        doc.connector(
                            ens_idx,
                            yp + 0.15,
                            cls=['output'])
                    )
                    show = True
    
                if show:
                    group.add(
                        doc.connector(
                            ens_idx,
                            yp,
                            txt,
                            cls=['unknown']
                        )
                    )

        group.translate(50, 0)

        doc.add(group)

        doc['class'] = 'movetree'

        left_x = -max_level * doc.scale_x - 120
        top_y = - 120
        width = len(self.ensembles) * doc.scale_x - left_x + 50
        height = (total + 1) * doc.scale_y - top_y

        # adjust viewbox to fit full image
        doc['viewBox'] = '%.2f %.2f %.2f %.2f' % (
            left_x,
            top_y,
            width,
            height
        )
        doc['width'] = width

        return doc


class PathTreeBuilder(Builder):
    def __init__(self):

        super(PathTreeBuilder, self).__init__(['movers'])
        self.rejected = False
        self.obj = list()
        self.doc = None

        self.move_list = {}
        self.step_list = {}
        self.samp_list = {}

        self.css_style = vis_css

        self.states = {}
        self.op = None

        self._sample_list = None
        self._steps = None

        self.options.movers.update({
            paths.ReplicaExchangeMover: {
                'name': 'RepEx',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'RepEx',
                'suffix': 'x',
                'label_position': 'left',
                'cls': ['repex']
            },
            paths.BackwardShootMover: {
                'name': 'Shooting',
                'overlap': 'none',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': '',
                'suffix': 'b',
                'label_position': 'left',
                'cls': ['shooting']
            },
            paths.ForwardShootMover: {
                'name': 'Shooting',
                'overlap': 'none',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': '',
                'suffix': 'f',
                'label_position': 'right',
                'cls': ['shooting']
            },
            paths.BackwardExtendMover: {
                'name': 'Extend',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'Extend',
                'suffix': 'b',
                'label_position': 'left',
                'cls': ['extend']
            },
            paths.ForwardExtendMover: {
                'name': 'Extend',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'Extend',
                'suffix': 'f',
                'label_position': 'right',
                'cls': ['extend']
            },
            paths.FinalSubtrajectorySelectMover: {
                'name': 'Truncate',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'Trunc',
                'suffix': 't',
                'label_position': 'right',
                'cls': ['extend']
            },
            paths.FirstSubtrajectorySelectMover: {
                'name': 'Truncate',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'Trunc',
                'suffix': 't',
                'label_position': 'left',
                'cls': ['extend']
            },
            paths.EnsembleHopMover: {
                'name': 'hop',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': 'EnsembleHop',
                'suffix': 'h',
                'label_position': 'left',
                'cls': ['hop']
            },
            paths.PathReversalMover: {
                'name': 'reversal',
                'overlap': 'line',
                'fw': '',
                'bw': '',
                'overlap_label': 'Reversal',
                'suffix': 'r',
                'label_position': 'left',
                'cls': ['reversal']
            },
            'new': {
                'name': 'new',
                'overlap': 'blocks',
                'fw': 'blocks',
                'bw': 'blocks',
                'suffix': '+',
                'overlap_label': 'New',
                'label_position': 'left',
                'cls': ['unknown']
            },
            'unknown': {
                'name': '???',
                'overlap': 'line',
                'fw': 'blocks',
                'bw': 'blocks',
                'overlap_label': '???',
                'suffix': '?',
                'label_position': 'left',
                'cls': ['repex']
            }
        })
        self.options.ui.update({
            'step': True,
            'correlation': True,
            'sample': True,
            'virtual': False,
            'cv': True,
            'info': False
        })
        self.options.analysis.update({
            'time_symmetric': True,
            'flip_time_direction': False,
            'joined_blocks': False
        })
        self.options.css.update({
            'scale_x': 5,
            'scale_y': 10,
            'zoom': 1.0,
            'horizontal_gap': False,
            'width': 'inherit'
        })
        self.options.format.update({
            'default_label': lambda x: hex(id(x))[-5:] + ' ',
            'trajectory_label': None,
            'sample_label': None,
            'step_label': None,
            'snapshot_label': None,
            'display_repeated': True,
            'new_snapshots': True,
            'repeated_snapshots': True
        })

    @staticmethod
    def construct_heritage(sample):
        return list(reversed(list(sample.heritage)))

    @property
    def samples(self):
        return list(self._sample_list)

    @samples.setter
    def samples(self, samples):
        self._sample_list = SampleList(samples)

    @property
    def steps(self):
        return self._steps

    @steps.setter
    def steps(self, steps):
        self._steps = steps

    def render(self):
        samples = self.samples
        doc = TreeRenderer(self.css_style)
        self.doc = doc

        opts = self.options

        doc.scale_x = opts.css['scale_x']
        doc.scale_y = opts.css['scale_y']
        if type(opts.css['horizontal_gap']) is bool:
            doc.horizontal_gap = 0.05 if opts.css['horizontal_gap'] else 0.0
        else:
            doc.horizontal_gap = opts.css['horizontal_gap']

        assume_reversed_as_same = opts.analysis['time_symmetric']

        trj_format = opts.format['trajectory_label'] or opts.format['default_label'] or (lambda obj: '')
        smp_format = opts.format['sample_label'] or opts.format['default_label'] or (lambda obj: '')
        snp_format = opts.format['snapshot_label'] or opts.format['default_label'] or (lambda obj: '')

        if opts.ui['info']:
            doc.defs.add(doc.script(
                content='''
                   box = $('.opstree .infobox text')[0];
                   var kernel = IPython.notebook.kernel;
                   $('.opstree .block').each(
                    function() {
                    json = JSON.parse($(this)[0].firstChild.textContent);
                     $(this).data(json);
                    }
                   );
                   $('.opstree .block').hover(
                    function(){
                      box.textContent =
                      'Snapshot(' + $(this).data('snp') + ')' + ' ' +
                      'Trajectoy(' + $(this).data('trj') + ')';
                     },
                     function(){
                      box.textContent = '';
                     });
            '''))

        # display_repeated = opts.format['display_repeated']
        # new_snapshots = opts.format['new_snapshots']
        # repeated_snapshots = opts.format['repeated_snapshots']

        group = doc.g(
            class_='tree'
        )

        matrix = self._sample_list.matrix
        
        for pos_y, sample in enumerate(samples):
            info = self.samp_list[sample]
            mover_type = info['mover_type']
            new_sample = info['new']
            shift = info['shift']
            length_fw = info['length_fw']
            length_bw = info['length_bw']
            length = info['length']
            # length_shared = info['length_shared']
            time_direction = info['time_direction']

            # overlap_reversed = info['overlap_reversed']
            # correlation = info['correlation']

            bw_cls = 'bw'
            fw_cls = 'fw'

            traj = sample.trajectory

            if time_direction == -1:
                traj = paths.Trajectory(list(reversed(list(traj))))
                bw_cls, fw_cls = fw_cls, bw_cls

            if new_sample:
                view_options = opts.movers['new']
            elif mover_type in opts.movers:
                view_options = opts.movers[mover_type]
            else:
                view_options = opts.movers['unknown']

            traj_str = str(trj_format(traj)) + view_options['suffix'].upper()

            cls = [] + view_options['cls']

            if sample in self.move_list:
                move = self.step_list[sample].change
                accepted = move.accepted
                if not accepted:
                    cls += ['rejected']

            if view_options['label_position'] == 'left':
                group.add(
                    doc.label(shift, pos_y, traj_str, cls=cls + ['left'])
                )
            elif view_options['label_position'] == 'right':
                group.add(
                    doc.label(shift + len(traj) - 1, pos_y, traj_str,
                              cls=cls + ['right'])
                )

            bw_x = shift + length_bw
            fw_x = shift + length - 1 - length_fw

            if 0 < length_bw:
                root_y = matrix.root(pos_y, bw_x)

                if root_y < pos_y:
                    group.add(
                        doc.vertical_connector(bw_x, root_y, pos_y,
                                               cls=cls + [bw_cls, 'connection'])
                    )

            if 0 < length_fw:
                root_y = matrix.root(pos_y, fw_x)

                if root_y < pos_y:
                    group.add(
                        doc.vertical_connector(fw_x, root_y, pos_y,
                                               cls=cls + [fw_cls, 'connection'])
                    )

            for add_cls, region, label in [
                    [[bw_cls], (0, length_bw), ''],
                    [[fw_cls], (length - length_fw, length), ''],
                    [['overlap', 'whiteback'], (length_bw, length - length_fw), view_options['overlap_label']],
                    [['reversed'], (length_bw, length - length_fw)]
            ]:
                vis_type = ''
                if vis_type == 'line':
                    group.add(
                        doc.horizontal_region(shift + region[0], pos_y, region[1] - region[0],
                                              view_options['overlap_label'], cls=cls + add_cls)
                    )
                elif vis_type == 'block':
                    group.add(
                        doc.block(
                            shift + region[0],
                            pos_y,
                            view_options['overlap_label'],
                            w=region[1] - region[0],
                            extend_left=False,
                            cls=cls + add_cls
                        ))
                elif vis_type == 'single':
                    for pos in range(region[0], region[1]):
                        snapshot = traj[pos]
                        pos_x = shift + pos

                        if opts.ui['info']:
                            data = {
                                'smp': smp_format(sample),
                                'snp': snp_format(snapshot),
                                'trj': trj_format(sample.trajectory)
                            }
                        else:
                            data = {}

                        txt = ''

                        if self.op is not None and opts.ui['cv']:
                            txt = str(self.op(snapshot))

                        group.add(
                            doc.block(
                                pos_x,
                                pos_y,
                                txt,
                                extend_left=pos > 0,
                                extend_right=pos < len(traj) - 1,
                                cls=cls + add_cls,
                                data=data
                            ))

        min_x, max_x = min(matrix.matrix.keys()), max(matrix.matrix.keys())
        min_y, max_y = 0, len(samples) - 1

        # print min_x, max_x
        # print min_y, max_y
        # print len(matrix), len(matrix[0])

        if hasattr(self, 'states') and self.states:
            for color, op in self.states.iteritems():
                xp = None
                for y in range(0, max_y - min_y + 1):
                    left = None
                    yp = y + min_y
                    for x in range(0, (max_x - min_x + 1)):
                        xp = x + min_x

                        if matrix[y, x] is not None and bool(op(matrix[y, x])):
                            if left is None:
                                left = xp
                        else:
                            if left is not None:
                                group.add(
                                    doc.shade(left, yp, xp - left, color=color)
                                )
                                left = None

                    if left is not None:
                        group.add(
                            doc.shade(left, yp, xp - left + 1, color=color)
                        )

        group.translate(32 + doc.w(1 - min_x), doc.h(1))

        tree_group = group

        group = doc.g(
            class_='legend'
        )

        if opts.ui['info']:
            group.add(
                doc.label(0, -1, 'Information', cls=['infobox'])
            )

        columns = 0
        tree_scale = opts.css['scale_x']
        doc.scale_x = 32

        if opts.ui['correlation']:
            columns += 1
            cor_x = -columns
        else:
            cor_x = None

        if opts.ui['sample']:
            columns += 1
            smp_x = -columns
        else:
            smp_x = None

        if opts.ui['step']:
            columns += 1
            cyc_x = -columns
        else:
            cyc_x = None

        if smp_x is not None:
            group.add(
                doc.label(smp_x, 0, 'smp')
            )

        if cyc_x is not None:
            group.add(
                doc.label(cyc_x, 0, 'cyc')
            )

        if cor_x is not None:
            group.add(
                doc.label(cor_x, 0, 'cor')
            )

        old_tc = 1

        width = 64 + tree_scale * (max_x - min_x + 2) - doc.scale_x * (-0.5 - columns)
        height = doc.scale_y * (max_y + 3.0)
        left_x = (-0.5 - columns) * doc.scale_x
        top_y = -1.5 * doc.scale_y

        if len(samples) > 0:
            prev = samples[0].trajectory
            cls = ['tableline']

            for tc, s in enumerate(samples):
                group.add(
                    doc.rect(
                        class_=doc.c(cls),
                        insert=doc.xy(-0.5 - columns, 1 + tc - 0.45),
                        size=(
                            width,
                            doc.scale_y * 0.9
                        )
                    )
                )
                if tc > 0:
                    if not paths.Trajectory.is_correlated(s.trajectory, prev, time_reversal=assume_reversed_as_same):
                        if cor_x is not None:
                            group.add(
                                doc.vertical_region(
                                    cor_x,
                                    old_tc,
                                    1 + tc - old_tc,
                                    cls=['correlation']
                                )
                            )

                        old_tc = 1 + tc
                        prev = s.trajectory

                if smp_x is not None:
                    group.add(
                        doc.label(smp_x, 1 + tc, str(
                            trj_format(s)))
                    )

                if cyc_x is not None:
                    if s in self.step_list:
                        txt = str(self.step_list[s].mccycle)
                    else:
                        txt = '---'

                    group.add(
                        doc.label(cyc_x, 1 + tc, str(
                            txt))
                    )

        if cor_x is not None:
            group.add(
                doc.vertical_region(
                    cor_x,
                    old_tc,
                    1 + len(samples) - old_tc,
                    extend_bottom=False,
                    cls=['correlation']))

        group_all = doc.g()
        group_all.add(group)
        group_all.add(tree_group)

        zoom = opts.css['zoom']

        group_all.scale(zoom)

        doc.add(group_all)

        # set the overall OPS tree class
        doc['class'] = 'opstree'

        # adjust viewbox to fit full image
        doc['viewBox'] = '%.2f %.2f %.2f %.2f' % (
            left_x * zoom,
            top_y * zoom,
            width * zoom,
            height * zoom
        )

        # set width
        w_opt = opts.css['width']
        if w_opt == 'inherit':
            doc['width'] = width * zoom
        else:
            doc['width'] = w_opt

        return doc

    def use_storage_indices(self, storage):
        self.options.format['default_label'] = storage.idx


class ReplicaHistoryTree(PathTreeBuilder):
    """
    Simplified PathTreeBuilder for the common case of tracking a replica
    over some steps.

    Intended behaviors:
    * The samples are determined during initialization.
    * The defaults are as similar to the old tree representation as
      reasonable.
    * This object also calculates decorrelated trajectories (which is
      usually what we look for from this tree). The number of decorrelated
      trajectories is obtained as the length of that list, and does not
      require an extra method.
    """

    def __init__(self, steps, replica):
        # steps variable here and just iterate over storage.
        super(ReplicaHistoryTree, self).__init__()
        self.replica = replica
        self._accepted_samples = None
        self._trial_samples = None

        # defaults:
        self.rejected = False
        self.states = []

        self.steps = steps

    @property
    def accepted_samples(self):
        """
        Returns the accepted samples in self.steps involving self.replica
        """
        if self._accepted_samples is None:
            samp = self.steps[-1].active[self.replica]
            samples = [samp]
            while samp.parent is not None:
                samp = samp.parent
                samples.append(samp)

            self._accepted_samples = list(reversed(samples))

        return self._accepted_samples

    @property
    def trial_samples(self):
        """
        Returns trial samples from self.steps involving self.replica
        """
        if self._trial_samples is None:
            samp = self.steps[0].active[self.replica]
            samples = [samp]
            for step in self.steps:
                rep_trials = [s for s in step.change.trials
                              if s.replica == self.replica]
                if len(rep_trials) > 0:
                    samples.append(rep_trials[-1])

            self._trial_samples = samples

        return self._trial_samples

    @property
    def decorrelated_trajectories(self):
        """List of decorrelated trajectories from the internal samples.

        In path sampling, two trajectories are said to be "decorrelated" if
        they share no frames in common. This is particularly important in
        one-way shooting. This function returns the list of trajectories,
        making the number (i.e., the length of the list) also easily
        accessible.
        """
        return None


class SnapshotMatrix(object):
    def __init__(self, sample_list):
        self.sample_list = sample_list
        self.matrix = {}
        self.shift = [0] * range(len(sample_list))

    def __setitem__(self, key, value):
        y_pos = key[0]
        x_pos = key[1]

        if x_pos not in self.matrix:
            self.matrix[x_pos] = {}

        if isinstance(value, paths.BaseSnapshot):
                self.matrix[x_pos][y_pos] = value

        elif type(value) is paths.Trajectory:
            for pos, snapshot in enumerate(value):
                self[y_pos, x_pos + pos] = snapshot

            self.shift[y_pos] = x_pos

    def __getitem__(self, item):
        y_pos = item[0]
        x_pos = item[1]
        if x_pos in self.matrix:
            return self.matrix[x_pos][y_pos]
        else:
            raise KeyError('Position not found.')

    def get(self, y_pos, x_pos):
        if x_pos in self.matrix:
            return self.matrix[x_pos].get(y_pos)
        else:
            return None

    def is_new(self, y_pos, x_pos):
        snapshot = self[y_pos, x_pos]

        x = self.matrix[x_pos]

        pos = y_pos
        while pos > 0:
            new_y_pos = self.sample_list.parent(pos)

            if not new_y_pos or new_y_pos > pos:
                return True

            pos = new_y_pos

            if snapshot is x[pos]:
                return False

        return True

    def root(self, y_pos, x_pos):
        snapshot = self[y_pos, x_pos]

        x = self.matrix[x_pos]

        pos = y_pos
        while pos > 0:
            new_y_pos = self.sample_list.parent(pos)

            if not new_y_pos or new_y_pos > pos:
                return pos

            if snapshot is not x[new_y_pos]:
                return pos

            pos = new_y_pos

        return pos

    def parent(self, y_pos, x_pos):
        snapshot = self[y_pos, x_pos]

        x = self.matrix[x_pos]

        if y_pos == 0:
            return None

        new_y_pos = self.sample_list.parent(y_pos)

        if not new_y_pos or new_y_pos > y_pos:
            return None

        if snapshot is not x[new_y_pos]:
            return

        return new_y_pos


class SampleList(list):
    """
    A timely ordered series of `Sample` objects.

    This is effectively a list object enhanced with a few additional functions that
    simplify analysis. Although this can hold an arbitrary list of samples it is meant
    to represent a time evolution of samples and thus samples that have a causal relation.

    Examples would be the history of samples that lead to a specific samples (heritage)
    or the history of samples in a specific ensemble or of a given replica.

    Last it provides some useful filters that make sense for samples.
    """

    def __init__(self, samples):
        list.__init__(self, samples)

        self.time_symmetric = True
        self.flip_time_direction = False
        self.matrix = []
        self.parents = []
        self.samp_list = {}

    def __add__(self, other):
        return SampleList(list.__add__(self, other))

    def __getitem__(self, item):
        if type(item) is slice:
            return SampleList(list.__getitem__(self, item))
        elif hasattr(item, '__iter__'):
            SampleList(list.__getitem__(self, item))
        else:
            return list.__getitem__(self, item)

    def parent(self, idx):
        try:
            if type(idx) is int:
                return self.index(self[idx].parent)
            else:
                return self.index(idx.parent)
        except ValueError:
            return None

    def _trajectory_index(self, trajectory, snapshot):
        if self.time_symmetric:
            return trajectory.index_symmetric(snapshot)
        else:
            return trajectory.index(snapshot)

    def _trajectory_contains(self, trajectory, snapshot):
        if self.time_symmetric:
            return trajectory.contains_symmetric(snapshot)
        else:
            return snapshot in trajectory

    def snapshot_position_x(self, sample, snapshot):
        if type(sample) is int:
            sample = self[sample]

        if sample in self.samp_list:
            if self.samp_list[sample]['time_direction'] > 0:
                x_pos = self.samp_list[sample]['shift'] + self._trajectory_index(sample, snapshot)
            else:
                x_pos = self.samp_list[sample]['shift'] + len(sample) - 1 - self._trajectory_index(sample, snapshot)

            return x_pos
        else:
            return None

    def analyze(self):

        matrix = SnapshotMatrix(self)
        samples = self

        flip_time_direction = self.flip_time_direction

        parent = None

        for y_pos, sample in enumerate(samples):
            mover_type = type(sample.mover)
            traj = sample.trajectory
            parent_shift = 0
            parent_traj = None
            overlap = None

            if sample.parent is not None:
                parent = sample.parent

            if parent not in self.samp_list:
                parent = None
                time_direction = +1

            if parent is not None:
                parent_shift = self.samp_list[parent]['shift']
                time_direction = self.samp_list[parent]['time_direction']

                parent_traj = parent.trajectory

                if time_direction == -1:
                    traj = paths.Trajectory(list(reversed(list(traj))))
                    parent_traj = paths.Trajectory(list(reversed(list(parent_traj))))

                overlap = parent_traj.shared_subtrajectory(traj, time_reversal=self.time_symmetric)
                overlap_length = len(overlap)

            if overlap is None or len(overlap) == 0:
                # no overlap so we need to start new
                traj_shift = 0

                self.samp_list[sample] = {
                    'shift': 0,
                    'new': True,
                    'mover_type': mover_type,
                    'time_direction': time_direction,
                    'correlation': 0.0
                }

            else:
                new_fw = self._trajectory_index(traj, overlap[-1])
                new_bw = self._trajectory_index(traj, overlap[0])

                overlap_reversed = False

                if new_bw > new_fw:
                    overlap_reversed = True

                    new_fw, new_bw = new_bw, new_fw

                    if flip_time_direction:
                        # reverse the time and adjust the shifting

                        traj = paths.Trajectory(list(reversed(list(traj))))
                        time_direction *= -1
                    else:
                        # after
                        overlap_length = 0

                traj_shift = parent_shift + self._trajectory_index(parent_traj, overlap[0]) - new_bw

                self.samp_list[sample] = {
                    'shift': traj_shift,
                    'length_fw': len(traj) - 1 - new_fw,
                    'length_bw': new_bw,
                    'length_shared': overlap_length,
                    'length': len(traj),
                    'overlap_reversed': overlap_reversed,
                    'new': False,
                    'mover_type': mover_type,
                    'time_direction': time_direction,
                    'correlation': (1.0 * overlap_length) / len(traj),
                    'parent_y': self.parent(sample)
                }

            matrix[y_pos, traj_shift] = traj

            parent = sample

        self.matrix = matrix

    @property
    def correlation(self):
        return [self.samp_list[s]['correlation'] for s in self]

    @property
    def decorrelated_trajectories(self):
        """List of decorrelated trajectories from the internal samples.

        In path sampling, two trajectories are said to be "decorrelated" if
        they share no frames in common. This is particularly important in
        one-way shooting. This function returns the list of trajectories,
        making the number (i.e., the length of the list) also easily
        accessible.
        """
        prev = self[0].trajectory
        decorrelated = [prev]

        for s in [samp for samp in self]:
            if not s.trajectory.is_correlated(prev, self.time_symmetric):
                decorrelated.append(s.trajectory)
                prev = s.trajectory

        return decorrelated


vis_css = r"""
.opstree text, .movetree text {
    alignment-baseline: central;
    font-size: 10px;
    text-anchor: middle;
    font-family: Futura-CondensedMedium;
    font-weight: lighter;
    stroke: none !important;
}
.opstree text.shadow {
    stroke-width: 3;
    stroke: white !important;
}
.opstree text.bw.label {
    text-anchor: end;
}
.opstree text.fw.label {
    text-anchor: start;
}
.opstree .block text, .movetree .block text {
    fill: white !important;
    stroke: none !important;
}
.opstree g.block:hover rect {
    opacity: 0.5;
}
.opstree .repex {
    fill: blue;
    stroke: blue;
}
.opstree .extend {
    fill: blue;
    stroke: blue;
}
.opstree .truncate {
    fill: blue;
    stroke: blue;
}
.opstree .new {
    fill: black;
    stroke: black;
}
.opstree .unknown {
    fill: gray;
    stroke: gray;
}
.opstree .hop {
    fill: blue;
    stroke: blue;
}
.opstree .correlation {
    fill: black;
    stroke: black;
}
.opstree .shooting.bw {
    fill: green;
    stroke: green;
}
.opstree .shooting.fw {
    fill: red;
    stroke: red;
}
.opstree .reversal {
    fill: gold;
    stroke: gold;
}
.opstree .virtual {
    opacity: 0.1;
    fill:gray;
    stroke: none;
}
.opstree line {
    stroke-width: 2px;
}
.opstree .label {
    fill: black !important;
}
.opstree .h-connector {
/*                stroke-dasharray: 3 3; */
}
.opstree .rejected {
    opacity: 0.3;
}
.opstree .orange {
    fill: orange;
}
.tableline {
    fill: gray;
    opacity: 0.0;
}
.tableline:hover {
    opacity: 0.2;
}
.opstree .left.label .shift {
    transform: translateX(-24px);
}
.opstree .right.label .shift {
    transform: translateX(+24px);
}
.opstree .infobox text {
    text-anchor: start;
}

.movetree .label .shift {
    transform: translateX(-18px);
}

.movetree .label text {
    text-anchor: end;
}
.movetree .v-connector {
    stroke: black;
}
.movetree .v-hook {
    stroke: black;
}
.movetree .ensembles .head .shift {
    transform: translateY(0px) rotate(270deg) ;
}
.movetree .ensembles .head text {
    text-anchor: start;
}
.movetree .connector.input {
    fill: green;
}
.movetree .connector.output {
    fill: red;
}
.movetree .unknown {
    fill: gray;
}
"""

# css_file = os.path.join(os.path.dirname(__file__), 'vis.css')
#
# with open(css_file, 'r') as content_file:
#     vis_css = content_file.read()

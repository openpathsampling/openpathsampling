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
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()
        self.doc = None

        self.move_list = {}
        self.step_list = {}
        self.samp_list = {}

        self.css_style = vis_css

        self.states = []
        self.op = None

        self._samples = None
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
            'snapshot_label': None
        })

    @staticmethod
    def construct_heritage(sample):
        return list(reversed(list(sample.heritage)))

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, samples):
        self._samples = samples
        self.analyze()

    @property
    def steps(self):
        return self._steps

    @steps.setter
    def steps(self, steps):
        self._steps = steps

    def analyze(self):
        samples = self._samples
        self.move_list = {}
        self.step_list = {}
        self.matrix = {}

        opts = self.options
        assume_reversed_as_same = opts.analysis['time_symmetric']
        flip_time_direction = opts.analysis['flip_time_direction']

        if self.steps is not None:
            for step in self.steps:
                for ch in step.change:
                    if ch.samples is not None:
                        for trial in ch.samples:
                            self.step_list[trial] = step
                            self.move_list[trial] = ch

        self.samp_list = {}
        p_x = dict()

        time_direction = +1

        for sample in samples:
            mover_type = type(sample.mover)
            traj = sample.trajectory

            if time_direction == -1:
                traj = paths.Trajectory(list(reversed(list(traj))))

            overlap_reversed = False
            index_bw = None
            index_fw = None

            for snapshot in range(len(traj)):
                snap = traj[snapshot]
                if snap in p_x or assume_reversed_as_same and snap.reversed in p_x:
                    connect_bw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                    index_bw = snapshot
                    shift_bw = connect_bw - snapshot
                    break

            new_sample = False
            if index_bw is None:
                index_bw = 0
                index_fw = len(traj) - 1
                # no overlap, so skip
                new_sample = True
                shift_bw = 0
                shift_fw = 0
            else:
                for snapshot in range(len(traj) - 1, -1, -1):
                    snap = traj[snapshot]
                    if snap in p_x or (assume_reversed_as_same and snap.reversed in p_x):
                        connect_fw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                        index_fw = snapshot
                        shift_fw = connect_fw - snapshot
                        break

                if shift_bw != shift_fw:
                    overlap_reversed = True

                # now we know that the overlap is between (including) [connect_bw, connect_fw]
                # and the trajectory looks like [bw, ...] + [old, ...] + [fw, ...]
                # with bw
                # [0, ..., index_bw -1] + [index_bw, ..., index_fw] + [index_fw + 1, ..., len(traj) - 1]
                # both shift_fw and shift_bw always exist and are the same if the trajectory is extended
                # or truncated or shoot from. If the overlapping trajectory is reversed before extending
                # then we get
                # [0, ..., index_fw -1] + [index_fw, ..., index_bw] + [index_bw + 1, ..., len(traj) - 1]
                # this can be checked by index_bw > index_fw or shift_fw != shift_bw

            p_x = {}

            if flip_time_direction and overlap_reversed:
                # reverse the time and adjust the shifting

                index_bw = len(traj) - 1 - index_bw
                index_fw = len(traj) - 1 - index_fw

                shift = (shift_fw + shift_bw) / 2 + index_bw - (len(traj) - 1 - index_fw)
                traj = paths.Trajectory(list(reversed(list(traj))))

                time_direction *= -1

            else:
                shift = (shift_fw + shift_bw) / 2

            for pos, snapshot in enumerate(traj):
                pos_x = shift + pos
                p_x[snapshot] = pos_x

            self.samp_list[sample] = {
                'shift': shift,
                'index_fw': index_fw,
                'index_bw': index_bw,
                'overlap_reversed': overlap_reversed,
                'new_sample': new_sample,
                'mover_type': mover_type,
                'time_direction': time_direction
            }

    def _write_snapshot_block(self, traj, x, y, first, last, rejected):
        p_x = self._p_x
        p_y = self._p_y

        if rejected:
            return

        for pos, snapshot in enumerate(traj[first:last + 1]):
            if self.options.analysis['time_symmetric']:
                if (snapshot not in p_x and snapshot.reversed not in p_x) or \
                        snapshot in p_x and p_x[snapshot] != x + pos + first or \
                        snapshot.reversed in p_x and p_x[snapshot.reversed] != x + pos + first:
                    p_x[snapshot] = x + pos + first
                    p_y[snapshot] = y
            else:
                if snapshot not in p_x or p_x[snapshot] != x + pos + first:
                    p_x[snapshot] = x + pos + first
                    p_y[snapshot] = y

    def render(self):
        samples = self._samples
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
        join_blocks = opts.analysis['joined_blocks']

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

        p_x = dict()
        p_y = dict()

        self._p_x = p_x
        self._p_y = p_y

        min_range_x = 10000
        max_range_x = -10000

        group = doc.g(
            class_='tree'
        )
        
        for pos_y, sample in enumerate(samples):
            info = self.samp_list[sample]
            mover_type = info['mover_type']
            new_sample = info['new_sample']
            shift = info['shift']
            index_fw = info['index_fw']
            index_bw = info['index_bw']
            time_direction = info['time_direction']

            rejected = False

            bw_cls = 'bw'
            fw_cls = 'fw'

            if index_fw < index_bw:
                index_bw, index_fw = index_fw, index_bw

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
                    rejected = True

            if view_options['label_position'] == 'left':
                group.add(
                    doc.label(shift, pos_y, traj_str, cls=cls + ['left'])
                )
            elif view_options['label_position'] == 'right':
                group.add(
                    doc.label(shift + len(traj) - 1, pos_y, traj_str,
                              cls=cls + ['right'])
                )

            if 0 < index_bw < len(traj) - 1:
                if traj[index_bw] in p_y:
                    pos_y_old = p_y[traj[index_bw]]
                elif assume_reversed_as_same and traj[index_bw].reversed in p_y:
                    pos_y_old = p_y[traj[index_bw].reversed]
                else:
                    pos_y_old = None

                if pos_y_old is not None:
                    group.add(
                        doc.vertical_connector(shift + index_bw, pos_y_old, pos_y,
                                               cls=cls + [bw_cls, 'connection'])
                    )
            if 0 < index_fw < len(traj) - 1:
                if traj[index_fw] in p_y:
                    pos_y_old = p_y[traj[index_fw]]
                elif assume_reversed_as_same and traj[index_fw].reversed in p_y:
                    pos_y_old = p_y[traj[index_fw].reversed]
                else:
                    pos_y_old = None

                if pos_y_old is not None:
                    group.add(
                        doc.vertical_connector(shift + index_fw + 1, pos_y_old, pos_y,
                                               cls=cls + [fw_cls, 'connection'])
                    )

            if view_options['overlap'] in ['line'] or view_options['overlap'] in ['blocks'] and join_blocks:
                if view_options['overlap'] in ['line']:
                    group.add(
                        doc.horizontal_region(shift + index_bw, pos_y, index_fw - index_bw + 1,
                                              view_options['overlap_label'], cls=cls + ['overlap', 'whiteback'])
                    )
                else:
                    group.add(
                        doc.block(
                            shift + index_bw,
                            pos_y,
                            view_options['overlap_label'],
                            w=index_fw - index_bw + 1,
                            extend_left=False,
                            cls=cls + ['overlap', 'whiteback']
                        ))

                self._write_snapshot_block(traj, shift, pos_y, index_bw, index_fw, rejected)

            if view_options['bw'] in ['line'] or view_options['bw'] in ['blocks'] and join_blocks:
                if view_options['bw'] in ['line']:
                    group.add(
                        doc.horizontal_region(shift + 0, pos_y, index_bw, cls=cls + [bw_cls])
                    )
                else:
                    group.add(
                        doc.block(
                            shift + 0,
                            pos_y,
                            w=index_bw,
                            extend_left=False,
                            cls=cls + [bw_cls]
                        ))

                self._write_snapshot_block(traj, shift, pos_y, 0, index_bw - 1, rejected)

            if view_options['fw'] in ['line'] or view_options['fw'] in ['blocks'] and join_blocks:
                if view_options['fw'] in ['line']:
                    group.add(
                        doc.horizontal_region(
                            shift + index_fw + 1,
                            pos_y,
                            len(traj) - (index_fw + 1),
                            cls=cls + [fw_cls]
                        )
                    )
                else:
                    group.add(
                        doc.block(
                            shift + index_fw + 1,
                            pos_y,
                            w=len(traj) - (index_fw + 1),
                            extend_right=False,
                            cls=cls + [fw_cls]
                        ))

                self._write_snapshot_block(traj, shift, pos_y, index_fw + 1, len(traj) - 1, rejected)

            for pos, snapshot in enumerate(traj):
                pos_x = shift + pos

                if opts.ui['info']:
                    data = {
                        'smp': smp_format(sample),
                        'snp': snp_format(snapshot),
                        'trj': trj_format(sample.trajectory)
                    }
                else:
                    data = {}

                if snapshot not in p_y or True:
                    txt = ''

                    if self.op is not None and opts.ui['cv']:
                        txt = str(self.op(snapshot))

                    b_cls = []

                    if not join_blocks:
                        if view_options['bw'] == 'blocks' and pos < index_bw:
                            b_cls += [bw_cls]
                        elif view_options['fw'] == 'blocks' and pos > index_fw:
                            b_cls += [fw_cls]
                        elif view_options['overlap'] == 'blocks':
                            b_cls += ['overlap']

                        if len(b_cls) > 0:
                            group.add(
                                doc.block(
                                    pos_x,
                                    pos_y,
                                    txt,
                                    extend_left=pos > 0,
                                    extend_right=pos < len(traj) - 1,
                                    cls=cls + b_cls,
                                    data=data
                                ))
                            if not rejected:
                                p_x[snapshot] = pos_x
                                p_y[snapshot] = pos_y

                        else:
                            if opts.ui['virtual']:
                                group.add(
                                    doc.block(
                                        pos_x,
                                        pos_y,
                                        txt,
                                        extend_left=pos > 0,
                                        extend_right=pos < len(traj) - 1,
                                        cls=cls + ['virtual'],
                                        data=data
                                    ))

                if pos_x < min_range_x:
                    min_range_x = pos_x

                if pos_x > max_range_x:
                    max_range_x = pos_x

        self.p_x = p_x
        self.p_y = p_y

        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = 0, len(samples) - 1

        matrix = self._to_matrix(min_y=min_y, max_y=max_y)

        # print min_x, max_x
        # print min_y, max_y
        # print len(matrix), len(matrix[0])

        if hasattr(self, 'states') and len(self.states) > 0:
            for color, op in self.states.iteritems():
                xp = None
                for y in range(0, max_y - min_y + 1):
                    left = None
                    yp = y + min_y
                    for x in range(0, (max_x - min_x + 1)):
                        xp = x + min_x
                        val = None

                        if (matrix[y][x] is not None and bool(op(matrix[y][x]))):
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

        group.translate(32 + doc.w(1 - min_range_x), doc.h(1))

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

        if max_range_x == -10000:
            max_range_x = 0
            min_range_x = 0

        old_tc = 1

        width = 64 + tree_scale * (max_range_x - min_range_x + 2) - doc.scale_x * (-0.5 - columns)
        height = doc.scale_y * (max_y + 3.0)
        left_x = (-0.5 - columns) * doc.scale_x
        top_y = -1.5 * doc.scale_y

        if len(samples) > 0:
            prev = samples[0].trajectory
            cls = ['tableline']
            if sample in self.move_list:
                move = self.step_list[sample].change
                accepted = move.accepted
                if not accepted:
                    cls += ['rejected']

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

    @staticmethod
    def _get_min_max(d):
        if len(d) > 0:
            return min(d.values()), max(d.values())
        else:
            return 0, 0

    def _to_matrix(self, min_x=None, max_x=None, min_y=None, max_y=None):
        if min_x is None or max_x is None:
            min_x, max_x = self._get_min_max(self.p_x)

        if min_y is None or max_y is None:
            min_y, max_y = self._get_min_max(self.p_y)

        matrix = [[None] * (max_x - min_x + 1) for n in range(max_y - min_y + 1)]

        for s in self.p_x:
            px = self.p_x[s]
            py = self.p_y[s]
            matrix[py - min_y][px - min_x] = s

        return matrix

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

    def rebuild(self):
        self.analyze()

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
        prev = self.samples[0].trajectory
        decorrelated = [prev]
        # TODO: this should be restricted to accepted samples
        for s in [samp for samp in self.samples]:
            if not paths.Trajectory.is_correlated(s.trajectory, prev):
                decorrelated.append(s.trajectory)
                prev = s.trajectory

        return decorrelated


class SnapshotMatrix(object):
    def __init__(self, sample_list):
        self.sample_list = sample_list
        length = len(sample_list)
        self.matrix = [{} for n in range(length)]
        self.ranges = [{} for n in range(length)]
        self.shift = [0 for n in range(length)]

    def __setitem__(self, key, value):
        y_pos = key[0]
        x_pos = key[1]
        if type(value) is paths.Trajectory:
            for pos, snapshot in enumerate(value):
                if x_pos not in self.matrix[y_pos]:
                    self.matrix[y_pos][x_pos + pos] = snapshot
                else:
                    raise KeyError('Position already set. Can only be set once')

            self.shift[y_pos] = x_pos

        elif isinstance(value, paths.BaseSnapshot):
            if x_pos not in self.matrix[y_pos]:
                self.matrix[y_pos][x_pos] = value
            else:
                raise KeyError('Position already set. Can only be set once')

    def __getitem__(self, item):
        y_pos = item[0]
        x_pos = item[1]
        return self.matrix[y_pos][x_pos]

    def get(self, y_pos, x_pos):
        if x_pos in self.matrix[y_pos]:
            return self.matrix[y_pos]
        else:
            return None

    def is_new(self, y_pos, x_pos):
        if x_pos not in self.matrix[y_pos]:
            raise KeyError('No snapshot at this position')

        snapshot = self.matrix[y_pos][x_pos]

        pos = y_pos
        while pos > 0:
            if snapshot in self.matrix[pos - 1]:
                return False

            pos -= 1

        return True

    def first(self, y_pos, x_pos):
        snapshot = self.matrix[y_pos][x_pos]

        pos = y_pos
        if pos > 0:
            if self.sample_list[pos].parent is None:
                return None

            new_y_pos = self.sample_list.parent(pos)

            if new_y_pos > pos:
                return None

            if self.matrix[new_y_pos][x_pos] is snapshot:
                return None

            return new_y_pos
        else:
            return None

    def parent(self, y_pos, x_pos):
        snapshot = self.matrix[y_pos][x_pos]

        pos = y_pos
        found = None
        while pos > 0:
            if self.sample_list[pos].parent is None:
                return found

            new_y_pos = self.sample_list.parent(pos)

            if new_y_pos > pos:
                return found

            if self.matrix[new_y_pos][x_pos] is snapshot:
                return found

            found = new_y_pos
            pos = new_y_pos
            snapshot = snapshot.parent

        return found

    def all_new(self, y_pos):


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
        if type(idx) is int:
            return self.index(self[idx].parent)
        else:
            return self.index(idx.parent)

    def analyze(self):

        matrix = SnapshotMatrix(self)
        samples = self

        # erase all steps from first onward
        if first > 0:
            self.matrix = self.matrix[:first]

        for pos in range(len(samples) - first):
            self.matrix.append({})

        assume_reversed_as_same = self.time_symmetric
        flip_time_direction = self.flip_time_direction

        time_direction = +1

        for y_pos, sample in enumerate(samples[first:], first):
            mover_type = type(sample.mover)
            traj = sample.trajectory

            if time_direction == -1:
                traj = paths.Trajectory(list(reversed(list(traj))))

            overlap_reversed = False
            index_bw = None
            index_fw = None

            for snapshot in range(len(traj)):
                snap = traj[snapshot]
                if snap in p_x or assume_reversed_as_same and snap.reversed in p_x:
                    connect_bw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                    index_bw = snapshot
                    shift_bw = connect_bw - snapshot
                    break

            new_sample = False
            if index_bw is None:
                index_bw = 0
                index_fw = len(traj) - 1
                # no overlap, so skip
                new_sample = True
                shift_bw = 0
                shift_fw = 0
            else:
                for snapshot in range(len(traj) - 1, -1, -1):
                    snap = traj[snapshot]
                    if snap in p_x or (assume_reversed_as_same and snap.reversed in p_x):
                        connect_fw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                        index_fw = snapshot
                        shift_fw = connect_fw - snapshot
                        break

                if shift_bw != shift_fw:
                    overlap_reversed = True

                # now we know that the overlap is between (including) [connect_bw, connect_fw]
                # and the trajectory looks like [bw, ...] + [old, ...] + [fw, ...]
                # with bw
                # [0, ..., index_bw -1] + [index_bw, ..., index_fw] + [index_fw + 1, ..., len(traj) - 1]
                # both shift_fw and shift_bw always exist and are the same if the trajectory is extended
                # or truncated or shoot from. If the overlapping trajectory is reversed before extending
                # then we get
                # [0, ..., index_fw -1] + [index_fw, ..., index_bw] + [index_bw + 1, ..., len(traj) - 1]
                # this can be checked by index_bw > index_fw or shift_fw != shift_bw

            p_x = {}

            if flip_time_direction and overlap_reversed:
                # reverse the time and adjust the shifting

                index_bw = len(traj) - 1 - index_bw
                index_fw = len(traj) - 1 - index_fw

                shift = (shift_fw + shift_bw) / 2 + index_bw - (len(traj) - 1 - index_fw)
                traj = paths.Trajectory(list(reversed(list(traj))))

                time_direction *= -1

            else:
                shift = (shift_fw + shift_bw) / 2

            for pos, snapshot in enumerate(traj):
                pos_x = shift + pos
                p_x[snapshot] = pos_x
                self.matrix[y_pos][pos_x] = snapshot

            self.samp_list[sample] = {
                'shift': shift,
                'index_fw': index_fw,
                'index_bw': index_bw,
                'overlap_reversed': overlap_reversed,
                'new_sample': new_sample,
                'mover_type': mover_type,
                'time_direction': time_direction
            }

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
            if not paths.Trajectory.is_correlated(s.trajectory, prev):
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

import svgwrite as svg
from svgwrite.container import Group
import os
import openpathsampling as paths
import networkx as nx
import json
import matplotlib.pyplot as plt
import StringIO
from networkx.readwrite import json_graph


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

    def _x(self, x):
        return self._w(x)

    def _y(self, y):
        return self._h(y)

    def _w(self, y):
        return self.scale_x * y

    def _h(self, y):
        return self.scale_y * y

    def _xy(self, x, y):
        return self._x(x), self._y(y)

    def _wh(self, w, h):
        return self._w(w), self._h(h)

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
            insert=self._xy(x - 0.5 + padding, y - 0.3),
            size=self._wh(1.0 * w - 2 * padding, 0.6),
        ))

        if extend_left:
            group.add(self.circle(
                center=self._xy(x - 0.5, y),
                r=self._w(padding)
            ))
        if extend_right:
            group.add(self.circle(
                center=(self._xy(x + w - 0.5, y)),
                r=self._w(padding)
            ))

        if extend_top:
            group.add(self.circle(
                center=self._xy(x, y - 0.3),
                r=self._w(padding)
            ))
        if extend_bottom:
            group.add(self.circle(
                center=(self._xy(x + w - 1.0, y + 0.3)),
                r=self._w(padding)
            ))

        group.add(self.text(
            text=str(text)[:4],
            insert=self._xy(x + (w - 1.0) / 2.0, y)
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
            start=self._xy(x - 0.5 + padding, y),
            end=self._xy(x - 0.5 + w - padding, y)
        ))

        if extend_left:
            group.add(self.circle(
                center=self._xy(x - 0.5, y),
                r=self._w(padding)
            ))
            group.add(self.line(
                start=self._xy(x - 0.5, y - 0.3),
                end=self._xy(x - 0.5, y + 0.3)
            ))

        if extend_right:
            group.add(self.circle(
                center=(self._xy(x + w - 0.5, y)),
                r=self._w(padding)
            ))
            group.add(self.line(
                start=self._xy(x + w - 0.5, y - 0.3),
                end=self._xy(x + w - 0.5, y + 0.3)
            ))

        group.add(self.text(
            text=str(text),
            insert=self._xy(x + (w - 1.0) / 2.0, y - 0.3)
        ))

        return group

    def vertical_region(self, x, y, w=1.0, text="", align="middle",
                        extend_top=True, extend_bottom=True, cls=None):

        if cls is None:
            cls = list()

        cls += ['v-region']

        padding = self.horizontal_gap
        width = 0.3
        gap = 0.0

        group = Group(
            class_=self.c(cls)
        )

        group.add(self.line(
            start=self._xy(x, y - 0.5 + gap),
            end=self._xy(x, y + w - 1 + 0.5 - gap)
        ))

        if extend_top:
            group.add(self.circle(
                center=self._xy(x, y - 0.5 + gap),
                r=self._w(padding)
            ))
            group.add(self.line(
                start=self._xy(x - width, y - 0.5 + gap),
                end=self._xy(x + width, y - 0.5 + gap)
            ))

        if extend_bottom:
            group.add(self.circle(
                center=(self._xy(x, y + (w - 1.0) + 0.5 - gap)),
                r=self._w(padding)
            ))
            group.add(self.line(
                start=self._xy(x - width, y + w - 1.0 + 0.5 - gap),
                end=self._xy(x + width, y + w - 1.0 + 0.5 - gap)
            ))

        group.add(self.text(
            text=str(text),
            insert=self._xy(x - width, y + (w - 1.0) / 2.0)
        ))

        return group

    def shade(self, x, y, w, cls=None):
        if cls is None:
            cls = list()

        cls += ['shade']

        return self.rect(
            class_=self.c(cls),
            insert=self._xy(x - 0.5, y + 0.35),
            size=self._wh(w, 0.1)
        )

    def vertical_connector(self, x, y1, y2, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-connector']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self._xy(x - 0.5, y1 + padding),
            end=self._xy(x - 0.5, y2 - padding)
        )

    def vertical_hook(self, x1, y1, x2, y2, cls=None):
        if cls is None:
            cls = list()

        cls += ['v-hook']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self._xy(x1, y1 + padding + 0.3),
            end=self._xy(x2, y2 - padding - 0.3)
        )

    def horizontal_connector(self, x1, x2, y, cls=None):
        if cls is None:
            cls = list()

        cls += ['h-connector']

        padding = self.horizontal_gap

        return self.line(
            class_=self.c(cls),
            start=self._xy(x1 + 0.5 + padding, y),
            end=self._xy(x2 - 0.5, y)
        )

    def label(self, x, y, text, cls=None):
        if cls is None:
            cls = list()

        cls += ['label']

        group = self.g(
            class_=self.c(cls)
        )

        group.translate(self._x(x), self._y(y))

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

    def vertical_label(self, x, y, w, text, cls=None):
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
            insert=self._xy(x, y),
            size=self._wh(w, h),
        )

    def to_svg(self):
        return self.tostring()

    def to_html(self, svg=None):
        if svg is None:
            svg = self.to_svg()

        html = '<!DOCTYPE html><html style="margin:0px; padding:0px; width:100%;">' + \
               svg + '<body style="margin:0px; padding:0px;"></body></html>'

        return html

    def _height(self):
        return self._h(self.height) + self.margin * 2

    def _width(self):
        return self._w(self.width) + self.margin * 2

    def write_html(self, file_name='tree.html'):
        with open(file_name, 'w') as f:
            f.write(self.to_html())

    def save_pdf(self, file_name='tree.pdf'):
        with open('tree_xxx.html', 'w') as f:
            f.write(self.to_html())

        bash_command = "open tree_xxx.html " + file_name

        os.system(bash_command)

    def clear(self):
        self.obj = []


class MoveTreeBuilder(object):
    def __init__(self, storage=None):
        self.rejected = False
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()
        self.storage = storage

        yp = 0
        self.traj_ens_x = dict()
        self.traj_ens_y = dict()

        self.traj_repl_x = dict()
        self.traj_repl_y = dict()

        self.ens_x = list()
        self.repl_x = list()

        # css_file = os.path.join(os.path.dirname(__file__), 'vis.css')
        #
        # with open(css_file, 'r') as content_file:
        #     self.css_style = content_file.read()

        self.css_style = vis_css


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

        tree = path.depth_pre_order(lambda this: this, only_canonical=True)

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

        for yp, (level, sub_mp) in enumerate(path.depth_pre_order(lambda this: this, only_canonical=True)):
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
                            '',
                            cls=['input']
                        )
                    )
    
                    show = True


                if out_ens is None or None in out_ens or ens in out_ens:
                    group.add(
                        doc.connector(
                            ens_idx,
                            yp + 0.15,
                            '',
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


class PathTreeBuilder(object):
    def __init__(self, storage, op=None, states=None):
        self.rejected = False
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()
        self.storage = storage
        self.doc = None

        self.move_list = {}
        self.step_list = {}
        self.samp_list = {}

        # css_file = os.path.join(os.path.dirname(__file__), 'vis.css')
        #
        # with open(css_file, 'r') as content_file:
        #     self.css_style = content_file.read()

        self.css_style = vis_css

        self.op = op
        self.show_redundant = False
        if states is None:
            states = []
        self.states = states

        self.options = {
            'mover': {
                paths.ReplicaExchangeMover: {
                    'name': 'RepEx',
                    'overlap': 'line',
                    'fw': 'blocks',
                    'bw': 'blocks',
                    'all': 'hidden',
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
                    'all': 'hidden',
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
                    'all': 'hidden',
                    'overlap_label': '',
                    'suffix': 'f',
                    'label_position': 'right',
                    'cls': ['shooting']
                },
                paths.EnsembleHopMover: {
                    'name': 'hop',
                    'overlap': 'line',
                    'fw': 'blocks',
                    'bw': 'blocks',
                    'all': 'hidden',
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
                    'all': 'blocks',
                    'overlap_label': 'Reversal',
                    'suffix': 'r',
                    'label_position': 'left',
                    'cls': ['reversal']
                },
                'new': {
                    'name': 'new',
                    'overlap': 'line',
                    'fw': 'blocks',
                    'bw': 'blocks',
                    'all': 'blocks',
                    'suffix': '+',
                    'overlap_label': '',
                    'label_position': 'left',
                    'cls': ['unknown']
                },
                'unknown': {
                    'name': '???',
                    'overlap': 'line',
                    'fw': 'blocks',
                    'bw': 'blocks',
                    'all': 'hidden',
                    'overlap_label': 'RepEx',
                    'suffix': '?',
                    'label_position': 'left',
                    'cls': ['repex']
                }
            },
            'ui': {
                'trajectory': True,
                'step': True,
                'correlation': True,
                'sample': True,
                'virtual': True,
                'cv': True,
            },
            'settings': {
                'register_rejected': False,
                'time_symmetric': True,
            },
            'geometry': {
                'scale_x': 15,
                'scale_y': 20,
                'horizontal_gap': True
            }
        }

    @staticmethod
    def construct_heritage(sample):
        list_of_samples = []

        samp = sample

        while samp.parent is not None:
            # just one sample so use this
            list_of_samples.append(samp)
            samp = samp.parent

        # reverse to get origin first
        return [samp for samp in reversed(list_of_samples)]

    def set_samples(self, samples):
        self._samples = samples
        self.analyze()

    def analyze(self):
        samples = self._samples
        self.move_list = {}
        self.step_list = {}

        options = self.options
        assume_reversed_as_same = options['settings']['time_symmetric']

        if len(samples) == 0:
            # no samples, nothing to do
            # TODO: Raise an exception or just ignore and don't output anything?
            return

        for step in self.storage.steps:
            for ch in step.change:
                if ch.samples is not None:
                    for trial in ch.samples:
                        self.step_list[trial] = step
                        self.move_list[trial] = ch

        self.samp_list = {}
        p_x = dict()

        for sample in samples:
            mover_type = type(sample.mover)
            traj = sample.trajectory

            overlap_reversed = False
            index_bw = None

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
                shift = 0
            else:
                for snapshot in range(len(traj) - 1, -1, -1):
                    snap = traj[snapshot]
                    if snap in p_x or (assume_reversed_as_same and snap.reversed in p_x):
                        connect_fw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                        index_fw = snapshot
                        shift_fw = connect_fw - snapshot
                        break

                # now we know that the overlap is between (including) [connect_bw, connect_fw]
                # and the trajectory looks like [bw, ...] + [old, ...] + [fw, ...]
                # with bw
                # [0, ..., index_bw -1] + [index_bw, ..., index_fw] + [index_fw + 1, ..., len(traj) - 1]
                # both shift_fw and shift_bw always exist and are the same if the trajectory is extended
                # or truncated or shoot from. If the overlapping trajectory is reversed before extending
                # then we get
                # [0, ..., index_fw -1] + [index_fw, ..., index_bw] + [index_bw + 1, ..., len(traj) - 1]
                # this can be checked by index_bw > index_fw or shift_fw != shift_bw

                shift = (shift_bw + shift_fw) / 2

                if index_bw > index_fw:
                    index_bw, index_fw = index_fw, index_bw
                    overlap_reversed = True

                    #            print sample.mover.__class__.__name__, 0, index_bw, index_fw, len(traj)-1

            for pos, snapshot in enumerate(sample.trajectory):
                pos_x = shift + pos
                p_x[snapshot] = pos_x


            self.samp_list[sample] = {
                'shift': shift,
                'index_fw': index_fw,
                'index_bw': index_bw,
                'overlap_reversed': overlap_reversed,
                'new_sample': new_sample,
                'mover_type': mover_type
            }

    def to_svg(self):
        return self.render().tostring()

    def render(self):
        samples = self._samples
        doc = TreeRenderer(self.css_style)
        self.doc = doc

        doc.scale_x = self.options['geometry']['scale_x']
        doc.scale_y = self.options['geometry']['scale_y']
        doc.horizontal_gap = 0.05 if self.options['geometry']['horizontal_gap'] else 0.0

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

        # /#                  kernel.execute('tv.frame = ' + $(this).data('snp')); /#

        p_x = dict()
        p_y = dict()

        min_range_x = 10000
        max_range_x = -10000

        t_count = 1

        group = doc.g(
            class_='tree'
        )

        options = self.options
        assume_reversed_as_same = options['settings']['time_symmetric']

        for sample in samples:
            info = self.samp_list[sample]
            mover_type = info['mover_type']
            new_sample = info['new_sample']
            shift = info['shift']
            index_fw = info['index_fw']
            index_bw = info['index_bw']
            traj = sample.trajectory

            if new_sample:
                view_options = options['mover']['new']
            elif mover_type in options['mover']:
                view_options = options['mover'][mover_type]
            else:
                view_options = options['mover']['unknown']

            # print shift, index_bw, index_fw

            traj_str = str(self.storage.idx(sample.trajectory)) + view_options['suffix']

            cls = [] + view_options['cls']

            if sample in self.move_list:
                move = self.step_list[sample].change
                accepted = move.accepted
                if not accepted:
                    cls += ['rejected']

            if view_options['label_position'] == 'left':
                group.add(
                    doc.label(shift, t_count, traj_str, cls=cls + ['left'])
                )
            elif view_options['label_position'] == 'right':
                group.add(
                    doc.label(shift + len(traj) - 1, t_count, traj_str,
                              cls=cls + ['right'])
                )

            if index_bw > 0:
                group.add(
                    doc.vertical_connector(shift + index_bw, p_y[traj[index_bw]], t_count,
                                           cls=cls + ['bw', 'connection'])
                )
            if index_fw < len(traj) - 1:
                group.add(
                    doc.vertical_connector(shift + index_fw + 1, p_y[traj[index_fw]], t_count,
                                           cls=cls + ['fw', 'connection'])
                )

            if view_options['overlap'] == 'line':
                group.add(
                    doc.horizontal_region(shift + index_bw, t_count, index_fw - index_bw + 1,
                                          view_options['overlap_label'], cls=cls + ['overlap'])
                )
                for pos, snapshot in enumerate(sample.trajectory[index_bw:index_fw + 1]):
                    if snapshot not in p_x:
                        p_x[snapshot] = shift + pos + index_bw
                        p_y[snapshot] = t_count

            if view_options['bw'] == 'line':
                group.add(
                    doc.horizontal_region(shift + 0, t_count, index_bw, cls=cls + ['bw'])
                )
                for pos, snapshot in enumerate(sample.trajectory[0:index_bw]):
                    p_x[snapshot] = shift + pos
                    p_y[snapshot] = t_count

            if view_options['fw'] == 'line':
                group.add(
                    doc.horizontal_region(shift + index_fw + 1, t_count, len(traj) - (index_fw + 1), cls=cls + ['fw'])
                )
                for pos, snapshot in enumerate(sample.trajectory[index_fw + 1:]):
                    p_x[snapshot] = shift + pos + index_fw + 1
                    p_y[snapshot] = t_count

            for pos, snapshot in enumerate(sample.trajectory):
                pos_x = shift + pos
                pos_y = t_count

                data = {
                    'smp': self.storage.idx(sample),
                    'snp': self.storage.idx(snapshot),
                    'trj': self.storage.idx(sample.trajectory)
                }

                if not snapshot in p_y or True:
                    txt = ''

                    if self.op is not None:
                        txt = self.op(snapshot)

                    b_cls = []

                    if view_options['bw'] == 'blocks' and pos < index_bw:
                        b_cls += ['bw']
                    elif view_options['fw'] == 'blocks' and pos > index_fw:
                        b_cls += ['fw']
                    elif view_options['overlap'] == 'blocks':
                        b_cls += ['overlap']

                    if len(b_cls) > 0:
                        group.add(
                            doc.block(
                                pos_x,
                                pos_y,
                                txt,
                                extend_left=pos > 0,
                                extend_right=pos < len(sample.trajectory) - 1,
                                cls=cls + b_cls,
                                data=data
                            ))

                        p_x[snapshot] = pos_x
                        p_y[snapshot] = pos_y

                    else:
                        if self.options['ui']['virtual']:
                            group.add(
                                doc.block(
                                    pos_x,
                                    pos_y,
                                    txt,
                                    extend_left=pos > 0,
                                    extend_right=pos < len(sample.trajectory) - 1,
                                    cls=cls + ['virtual'],
                                    data=data
                                ))

                if pos_x < min_range_x:
                    min_range_x = pos_x

                if pos_x > max_range_x:
                    max_range_x = pos_x

            t_count += 1

        self.p_x = p_x
        self.p_y = p_y

        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = self._get_min_max(self.p_y)

        matrix = self._to_matrix()

        if hasattr(self, 'states') and len(self.states) > 0:
            for color, op in self.states.iteritems():
                xp = None
                for y in range(0, max_y - min_y + 1):
                    left = None
                    yp = y + min_y
                    for x in range(0, (max_x - min_x + 1)):
                        xp = x + min_x

                        if matrix[y][x] is not None and bool(op(matrix[y][x])):
                            if left is None:
                                left = xp
                        else:
                            if left is not None:
                                group.add(
                                    doc.shade(left, yp, xp - left, cls=[color])
                                )
                                left = None

                    if left is not None:
                        group.add(
                            doc.shade(left, yp, xp - left + 1, cls=[color])
                        )

        group.translate(32 + doc._w(1 - min_range_x), 0)

        tree_group = group

        group = doc.g(
            class_='legend'
        )

        group.add(
            doc.label(0, 0, 'Information', cls=['infobox'])
        )


        columns = 0
        tree_scale = self.options['geometry']['scale_x']
        doc.scale_x = 32

        if self.options['ui']['correlation']:
            columns += 1
            cor_x = -columns
        else:
            cor_x = None

        if self.options['ui']['sample']:
            columns += 1
            smp_x = -columns
        else:
            smp_x = None

        if self.options['ui']['step']:
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

        prev = samples[0].trajectory
        old_tc = 1

        width = 64 + tree_scale * (max_range_x - min_range_x + 2) - doc.scale_x * (-0.5 - columns)
        height = doc.scale_y * (max_y + 1.0)
        left_x = (-0.5 - columns) * doc.scale_x
        top_y = -0.5 * doc.scale_y

        for tc, s in enumerate(samples):
            group.add(
                doc.rect(
                    class_=doc.c(['tableline']),
                    insert=doc._xy(-0.5 - columns, 1 + tc - 0.45),
                    size=(
                        width,
                        doc.scale_y * 0.9
                    )
                )
            )
            if tc > 0:
                if not paths.Trajectory.is_correlated(s.trajectory, prev):
                    if cor_x is not None:
                        group.add(
                            doc.vertical_region(
                                cor_x,
                                old_tc,
                                1 + tc - old_tc,
                                "",
                                cls=['correlation']
                            )
                        )

                    old_tc = 1 + tc
                    prev = s.trajectory

            if smp_x is not None:
                group.add(
                    doc.label(smp_x, 1 + tc, str(
                        self.storage.idx(s)))
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
                    "",
                    extend_bottom=False,
                    cls=['correlation']))

        doc.add(group)
        doc.add(tree_group)

        # set the overall OPS tree class
        doc['class'] = 'opstree'

        # adjust viewbox to fit full image
        doc['viewBox'] = '%.2f %.2f %.2f %.2f' % (
            left_x,
            top_y,
            width,
            height
        )

        return doc

    def _get_min_max(self, d):
        return min(d.values()), max(d.values())

    def _to_matrix(self):
        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = self._get_min_max(self.p_y)

        matrix = [[None] * (max_x - min_x + 1) for n in
                  range(max_y - min_y + 1)]

        for s in self.p_x:
            px = self.p_x[s]
            py = self.p_y[s]
            matrix[py - min_y][px - min_x] = s

        return matrix


class SVGDiGraph(nx.DiGraph):
    def _repr_svg_(self):
        plt.ioff()  # turn off interactive mode
        fig = plt.figure(figsize=(2, 2))
        ax = fig.add_subplot(111)
        nx.shell(self, ax=ax)
        output = StringIO.StringIO()
        fig.savefig(output, format='svg')
        plt.ion()  # turn on interactive mode
        return output.getvalue()


class MoveTreeNX(object):
    """Class to create a networkX based representation of a pathmover and change

    """

    def __init__(self, pathmover):
        self.pathmover = pathmover
        self._G = None

    @property
    def _enumeration(self):
        return enumerate(self.pathmover.depth_post_order(lambda this: this))

    @property
    def G(self):
        if self._G is None:
            G = nx.DiGraph()

            node_list = dict()

            for idx, data in self._enumeration:
                level, node = data
                node_list[node] = idx
                G.add_node(idx, name=node.name)

            for idx, data in self._enumeration:
                level, node = data
                subnodes = node.submovers
                for subnode in subnodes:
                    G.add_edge(node_list[node], node_list[subnode])

            self._G = G

        return self._G

    def draw(self):
        G = self.G
        pos = nx.spring_layout(G)

        for idx, data in self._enumeration:
            level, node = data
            nx.networkx_nodes(
                G,
                pos,
                nodelist=[idx],
                node_color='r',
                node_size=500,
                alpha=0.8
            )

        plt.axis('off')
        plt.show()

        return G

    @property
    def json_tree(self):
        data = json_graph.tree_data(self.G, len(self.G) - 1)
        return json.dumps(data)

    @property
    def json_node_link(self):
        data = json_graph.node_link_data(self.G)
        return json.dumps(data)

    def d3vis(self):
        return '''
        <style>

        .node circle {
          fill: #fff;
          stroke: steelblue;
          stroke-width: 1.5px;
        }

        .node {
          font: 10px sans-serif;
          stroke: black;
          stroke-width:0.35px;
        }

        .link {
          fill: none;
          stroke: #ccc;
          stroke-width: 1.5px;
        }

        </style>
        <div><svg id="d3-circ-tree-svg"></svg></div>
        ''' + '<script>var graph = ' + self.json_tree + ';</script>' + \
               '''
        <script src="http://d3js.org/d3.v3.min.js"></script>

        <script>

//        require.config({paths: {d3: "http://d3js.org/d3.v3.min"}});

//        require(["d3"], function(d3) {
        (function() {
            var diameter = 1200;
            var padding = 0;

            var tree = d3.layout.tree()
                .size([360, diameter / 2 - 120])
                .separation(function(a, b) { return (a.parent == b.parent ? 1 : 2) / a.depth; });

            var nodes = tree.nodes(graph),
                links = tree.links(nodes);

            var diagonal = d3.svg.diagonal.radial()
                .projection(function(d) { return [d.y, d.x / 360 * Math.PI]; });

            var vertical = d3.svg.diagonal.radial()
                .projection(function(d) { return [-d.x, -d.y / 360 * Math.PI]; });

            var svg = d3.select("#d3-circ-tree-svg")
                .attr("width", (diameter + padding))
                .attr("height", (diameter + padding) / 2)
              .append("g")
                .attr("transform", "translate(" + (diameter + padding) / 2 + "," + (0*(diameter + padding) / 2 + 50) + ")rotate(90)");

            var link = svg.selectAll(".link")
              .data(links)
              .enter().append("path")
              .attr("class", "link")
              .attr("d", diagonal);

            link.append("circle")
              .attr("r", 4.5);

            var node = svg.selectAll(".node")
              .data(nodes)
            .enter().append("g")
              .attr("class", "node")
              .attr("transform", function(d) { return "rotate(" + (d.x / 2 - 90) + ")translate(" + d.y + ")"; })

            node.append("circle")
              .attr("r", 4.5);

            node.append("text")
              .attr("dy", ".31em")
              .attr("text-anchor", function(d) { return d.x < 180 ? "start" : "end"; })
              .attr("transform", function(d) { return d.x < 180 ? "translate(8)" : "rotate(180)translate(-8)"; })
              .textual(function(d) { return d.name; });

            d3.select(self.frameElement).style("height", diameter + 50 + "px");
//        });
        })();

        </script>
        '''


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

    def __init__(self, storage, steps, replica):
        # TODO: if we implement substorages (see #330) we can remove the
        # steps variable here and just iterate over storage.
        super(ReplicaHistoryTree, self).__init__(storage)
        self.replica = replica
        self.steps = steps
        self._accepted_samples = None
        self._trial_samples = None

        # defaults:
        self.rejected = False
        self.show_redundant = False
        self.states = []

        # build the tree 
        self.set_samples(self.samples)

    def rebuild(self):
        """Rebuild the internal structures.

        It seems like some changes in the visualization require a complete
        rebuild. That's not ideal. If that can be changed, this function
        could be removed.

        """
        self.set_samples(self.samples)

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
    def samples(self):
        if self.rejected:
            return self.trial_samples
        else:
            return self.accepted_samples

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

vis_css = r"""
.opstree text, .movetree text {
    alignment-baseline: central;
    font-size: 10px;
    text-anchor: middle;
    font-family: Futura-CondensedMedium;
    font-weight: lighter;
    stroke: none !important;
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
import svgwrite
import os

import openpathsampling as paths
import networkx as nx

import json
import matplotlib.pyplot as plt
import StringIO

from networkx.readwrite import json_graph

# CSS attributes
# text_anchor=align,
# alignment_baseline='middle',
# font_family=self.font_family,
# fill=color

class TreeRenderer(object):
    def __init__(self):
        self.start_x = 0
        self.start_y = 0
        self.scale_x = 24
        self.scale_y = 24
        self.scale_th = 24
        self.document = None
        self.horizontal_gap = 0.05
        self.stroke_width = 0.1
        self.min_x = 10000
        self.min_y = 10000
        self.max_x = -10000
        self.max_y = -10000

        self.shift_x = 0
        self.shift_y = 0
        self.height = 0
        self.width = 0

        self.obj = list()
        self.margin = 0
        self.zoom = 1.0

    def __getattr__(self, item):
        if item in ['block', 'shade', 'v_connection', 'h_connection', 'label', 'rect',
                    'connector', 'v_hook', 'vertical_label', 'text', 'range', 'h_range']:
            # This will delay the execution of the draw commands until
            # we know where to draw
            return self._delay(object.__getattribute__(self, 'draw_' + item))
        else:
            return object.__getattribute__(self, item)

    def reset_pad(self):
        self.min_x = 10000
        self.min_y = 10000
        self.max_x = -10000
        self.max_y = -10000

    @staticmethod
    def _delay(func):
        def wrapper(*args, **kwargs):
            def fnc():
                return func(*args, **kwargs)

            return fnc

        return wrapper

    def add(self, obj):
        self.obj.append(obj)

    def pre(self, obj):
        self.obj.insert(0, obj)

    def _x(self, x):
        return self.start_x + self._w(x - self.shift_x)

    def _y(self, y):
        return self.start_y + self._h(y + self.shift_y)

    def _w(self, y):
        return self.scale_x * self.zoom * y

    def _h(self, y):
        return self.scale_y * self.zoom * y

    def _th(self, th):
        return self.scale_th * self.zoom * th

    def _xy(self, x, y):
        return self._x(x), self._y(y)

    def _wh(self, w, h):
        return self._w(w), self._h(h)

    def _xb(self, x, y):
        return self._x(x), self._y(y)

    def _pad(self, x, y, w, h):
        self.min_x = min(self.min_x, x, x + w)
        self.min_y = min(self.min_y, y, y + h)
        self.max_x = max(self.max_x, x + w, x)
        self.max_y = max(self.max_y, y + h, y)

    def draw_connector(self, x, y, text="", cls=None):

        if cls is None:
            cls=list()

        cls += ['connector']

        return self.draw_block(x, y, text, False, False, True, True, cls=cls)


    def draw_block(self, x, y, text="",
                   extend_right=True, extend_left=True,
                   extend_top=False, extend_bottom=False,
                   w=1.0, cls=None):

        if cls is None:
            cls=list()

        cls += ['block']

        document = self.document
        padding = self.horizontal_gap

        self._pad(x - 0.5 - padding, y - 0.3, 1.0 * w + 2* padding, 0.6)

        ret = list()

        ret.append(document.rect(
            class_= self.c(cls),
            insert=self._xy(x - 0.5 + padding, y - 0.3),
            size=self._wh(1.0 * w - 2 * padding, 0.6),
        ))
        if extend_left:
            ret.append(document.circle(
                class_= self.c(cls),
                center=self._xy(x - 0.5, y),
                r=self._w(padding)
            ))
        if extend_right:
            ret.append(document.circle(
                class_= self.c(cls),
                center=(self._xy(x + w - 0.5, y)),
                r=self._w(padding)
            ))

        ret.append(document.text(
            class_= self.c(cls),
            text=str(text)[:4],
            insert=self._xb(x + (w - 1.0) / 2.0, y)
        ))

        return ret

    @staticmethod
    def c(cls):
        return ' '.join(cls)

    def draw_range(self, x, y, w=1.0, text="",
                   extend_right=True, extend_left=True, cls=None):

        if cls is None:
            cls=list()

        cls += ['range']

        if w == 0:
            return []

        document = self.document
        padding = self.horizontal_gap

        self._pad(x - 0.5, y - 0.3, 1.0 * w, 0.6)

        ret = list()

        ret.append(document.rect(
            class_=self.c(cls),
            insert=self._xy(x - 0.5 + padding, y - 0.05),
            size=self._wh(1.0 * w - 2 * padding, 0.1)
        ))

        if extend_left:
            ret.append(document.circle(
                class_=self.c(cls),
                center=self._xy(x - 0.5, y),
                r=self._w(padding)
            ))
            ret.append(document.rect(
                class_=self.c(cls),
                insert=self._xy(x - 0.5 + padding, y - 0.3),
                size=self._wh(0.1, 0.6)
            ))

        if extend_right:
            ret.append(document.circle(
                class_=self.c(cls),
                center=(self._xy(x + w - 0.5, y)),
                r=self._w(padding)
            ))
            ret.append(document.rect(
                class_=self.c(cls),
                insert=self._xy(x + w - 0.6 - padding, y - 0.3),
                size=self._wh(0.1, 0.6)
            ))

        ret.append(document.text(
            class_=self.c(cls),
            text=str(text),
            insert=self._xb(x + (w - 1.0) / 2.0, y - 0.3)
        ))

        return ret

    def draw_h_range(self, x, y, w=1.0, color="blue", text="", align="middle",
                   extend_top=True, extend_bottom=True, cls=None):

        if cls is None:
            cls=list()

        cls += ['h-range']

        document = self.document
        padding = self.horizontal_gap

        self._pad(x - 0.5, y - 0.3, 1.0 * w, 0.6)

        ret = list()

        ret.append(document.rect(
            class_=self.c(cls),
            insert=self._xy(x, y - 0.3),
            size=self._wh(0.1, 1.0 * w - 0.4)
        ))

        if extend_top:
            ret.append(document.circle(
                class_=self.c(cls),
                center=self._xy(x, y - 0.3),
                r=self._w(padding)
            ))
            ret.append(document.rect(
                class_=self.c(cls),
                insert=self._xy(x - 0.5 + padding, y - 0.3),
                size=self._wh(1.0 - 2.0 * padding, 0.1)
            ))

        if extend_bottom:
            ret.append(document.circle(
                class_=self.c(cls),
                center=(self._xy(x, y + (w - 1.0) + 0.3)),
                r=self._w(padding)
            ))
            ret.append(document.rect(
                class_=self.c(cls),
                insert=self._xy(x - 0.5 + padding, y + (w - 1.0) + 0.3 - 0.1),
                size=self._wh(1.0 - 2.0 * padding, 0.1)
            ))

        ret.append(document.text(
            class_=self.c(cls),
            text=str(text),
            insert=self._xb(x - 0.3, y  + (w - 1.0) / 2.0)
        ))

        return ret

    def draw_text(self, x, y, text="", cls=None):

        if cls is None:
            cls=list()

        cls += ['text']

        document = self.document
        self._pad(x - 0.5, y, 1.0, 5.0)

        ret = list()

        ret.append(document.text(
            class_=self.c(cls),
            text=str(text)[:4],
            insert=self._xb(x, y),
            transform='rotate(90deg)'
        ))

        return ret

    def draw_shade(self, x, y, w, cls=None):
        if cls is None:
            cls=list()

        cls += ['shade']

        document = self.document
        self._pad(x - 0.5, y - 0.35, 1.0, 0.7)

        return [document.rect(
            class_=self.c(cls),
            insert=self._xy(x - 0.5, y + 0.35),
            size=self._wh(w, 0.1)
        )]

    def draw_shade_alt(self, x, y, w, cls=None):
        if cls is None:
            cls=list()

        cls += ['shade-alt']

        document = self.document

        self._pad(x - 0.5, y - 0.35, 1.0, 0.7)

        return [document.rect(
            class_=self.c(cls),
            insert=self._xy(x - 0.5, y - 0.35),
            size=self._wh(w, 0.7),
        )]

    def draw_v_connection(self, x, y1, y2, cls=None):
        if cls is None:
            cls=list()

        cls += ['v-connection']

        document = self.document
        stroke_width = self.stroke_width
        padding = self.horizontal_gap

        self._pad(x - 0.5 - stroke_width, y1,
                  2 * stroke_width, y2 - y1)

        return [document.line(
            class_=self.c(cls),
            start=self._xy(x - 0.5, y1 + padding),
            end=self._xy(x - 0.5, y2 - padding)
        )]

    def draw_v_hook(self, x1, y1, x2, y2, cls=None):
        if cls is None:
            cls=list()

        cls += ['v-hook']

        document = self.document
        stroke_width = self.stroke_width
        padding = self.horizontal_gap

        self._pad(x1 - stroke_width, y1,
                  2 * stroke_width + x2 - x1, y2 - y1)

        return [document.line(
            class_=self.c(cls),
            start=self._xy(x1, y1 + padding + 0.3),
            end=self._xy(x2, y2 - padding - 0.3)
        )]

    def draw_h_connection(self, x1, x2, y, cls=None):
        if cls is None:
            cls=list()

        cls += ['h-connection']

        document = self.document
        stroke_width = self.stroke_width
        padding = self.horizontal_gap

        self._pad(x1, y - stroke_width,
                  x2 - x1, 2 * stroke_width)

        return [document.line(
            class_=self.c(cls),
            start=self._xy(x1 + 0.5 + padding, y),
            end=self._xy(x2 - 0.5, y)
        )]

    def draw_label(self, x, y, w, text, cls=None):
        if cls is None:
            cls=list()

        cls += ['label']

        document = self.document

        self._pad(
            x-w, y-0.4, w, 0.8
        )

        return [document.text(
            class_=self.c(cls),
            text=str(text),
            insert=self._xb(x, y),
        )]

    def draw_vertical_label(self, x, y, w, text, cls=None):
        if cls is None:
            cls=list()

        cls += ['v-label']

        document = self.document

        self._pad(
            x, y,
            w, 0.6
        )

        return [document.text(
            text=str(text),
            insert=(0,0),
            transform='translate(' +
                ' '.join(map(str, self._xb(x, y))) +
                ')rotate(270)'
        )]

    def draw_rect(self, x, y, w, h, cls=None):
        if cls is None:
            cls=list()

        return [self.document.rect(
            class_=self.c(cls),
            insert=self._xy(x, y),
            size=self._wh(w, h),
        )]

    def to_document(self):

        self.start_x = self.margin
        self.start_y = self.margin

        doc = svgwrite.Drawing()
        self.document = doc

        doc.defs.add(doc.style(
            self.style
        ))
#        doc['width'] = str(self._width()) + 'px'
#        doc['height'] = str(self._height()) + 'px'
#        doc['width'] = '100%'
#        doc['height'] = '100%'

        doc['class'] = 'opstree'

        for obj in self.obj:
            parts = obj()
            map(doc.add, parts)

        # adjust viewbox to fit full image
        doc['viewBox'] = '%.2f %.2f %.2f %.2f' % (
            self._xy(self.min_x, self.min_y) + self._wh(self.max_x - self.min_x, self.max_y - self.min_y)
        )

        return doc

    def to_svg(self):
        doc = self.to_document()
        return doc.tostring()

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
        h = self._height()
        w = self._width()

        h_margin = 3.5
        v_margin = 3.5

        page_height = str(25.4 / 75.0 * h + v_margin) + 'mm'
        page_width = str(25.4 / 75.0 * w + h_margin) + 'mm'

        with open('tree_xxx.html', 'w') as f:
            f.write(self.to_html())

        bash_command = "wkhtmltopdf -l --page-width " + page_width + \
                      " --page-height " + page_height + \
                      " --disable-smart-shrinking " + \
                      "-B 1mm -L 1mm -R 1mm -T 1mm tree_xxx.html " + file_name

        os.system(bash_command)

    def clear(self):
        self.obj = []


class MoveTreeBuilder(object):
    def __init__(self, storage=None, op=None, states=None):
        self.rejected = False
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()
        self.storage = storage
        self.renderer = TreeRenderer()
        self.op = op
        if states is None:
            states = {}
        self.states = states

        self.t_count = 0
        self.traj_ens_x = dict()
        self.traj_ens_y = dict()

        self.traj_repl_x = dict()
        self.traj_repl_y = dict()

        self.ens_x = list()
        self.repl_x = list()

    def full(self, ensembles, clear=True):

        storage = self.storage

        if clear:
            self.renderer.clear()
        level_y = dict()

        old_sset = paths.SampleSet([])

        self.t_count = 1

        self.ens_x = [None] * len(ensembles)
        self.repl_x = [None] * len(ensembles)

        for sset in storage.samplesets[0:2]:
            path = sset.movepath
            # level_y = dict()

            for level, sub in path.depth_post_order(lambda this: tuple(
                    [this, old_sset.apply_samples(this.samples)])):
                self.t_count += 1

                sub_mp, sub_set = sub

                if sub_mp.__class__ is paths.SamplePathMoveChange:
                    self.renderer.add(
                        self.renderer.block(-8.0 + level, self.t_count, 'blue'))
                    self.renderer.add(
                        self.renderer.label(
                            -8.0 + level, self.t_count, 3,
                            sub_mp.mover.__class__.__name__[:-5],
                            align='start', color='black')
                    )
                else:
                    self.renderer.add(
                        self.renderer.block(-8.0 + level, self.t_count,
                                            'green'))
                    self.renderer.add(
                        self.renderer.label(-8.0 + level, self.t_count, 3,
                                            sub_mp.__class__.__name__[:-8],
                                            align='start', color='black')
                    )

                if level + 1 in level_y \
                        and level_y[level + 1] == self.t_count - 1:
                    self.renderer.add(
                        self.renderer.v_connection(-7.0 + level, self.t_count,
                                                   self.t_count - 1, 'black')
                    )
                    del level_y[level + 1]

                if level in level_y and level_y[level]:
                    self.renderer.add(
                        self.renderer.v_connection(-8.0 + level, self.t_count,
                                                   level_y[level], 'black')
                    )

                level_y[level] = self.t_count

                self.render_ensemble_line(ensembles, sub_set, color='gray')
                self.render_replica_line(len(ensembles), sub_set, color='gray')

            self.t_count += 1

            self.render_ensemble_line(ensembles, sset)
            self.render_replica_line(len(ensembles), sset)

            self.renderer.add(self.renderer.block(-8.0, self.t_count, 'black'))
            self.renderer.add(
                self.renderer.label(-8.0, self.t_count, 3,
                                    'storage.sampleset[%d]' % sset.idx[storage.samplesets],
                                    align='start', color='black')
            )

            old_sset = sset
            self.t_count += 1

        self.renderer.shift_x = - 9.0
        self.renderer.shift_y = -0.5
        self.renderer.height = 1.0 * self.t_count + 1.0
        self.renderer.width = 1.0 * len(ensembles) + 20.5


    def mover(self, pathmover, ensembles, clear=True):

        storage = self.storage

        if clear:
            self.renderer.clear()
        level_y = dict()

        self.t_count = 1

        self.ens_x = [None] * len(ensembles)
        self.repl_x = [None] * len(ensembles)

        path = pathmover

        self.render_ensemble_mover_head(ensembles, 5)

        self.t_count = 5

        total = len(path)

        for ens_idx, ens in enumerate(ensembles):
            self.renderer.add(
                self.renderer.v_hook( ens_idx, self.t_count, ens_idx,
                                           self.t_count + total + 1, 'gray')
            )

        for level, sub in path.depth_pre_order(lambda this: tuple(
                [this, None])):
            self.t_count += 1

            sub_mp, sub_set = sub

            name = sub_mp.name

            if sub_mp.__class__ is paths.SamplePathMoveChange:
                self.renderer.add(
                    self.renderer.block(-8.0 + level, self.t_count, 'blue'))
                self.renderer.add(
                    self.renderer.label(
                        -8.0 + level, self.t_count, 3,
                        sub_mp.mover.__class__.__name__[:-5],
                        align='start', color='black')
                )
            else:
                self.renderer.add(
                    self.renderer.block(-8.0 + level, self.t_count,
                                        'green'))
                self.renderer.add(
                    self.renderer.label(-8.0 + level, self.t_count, 3,
                                        sub_mp.__class__.__name__[:-5],
                                        align='start', color='black')
                )

            if False:
                if level + 1 in level_y \
                        and level_y[level + 1] == self.t_count - 1:
                    self.renderer.add(
                        self.renderer.v_connection(-7.0 + level, self.t_count,
                                                   self.t_count - 1, 'orange')
                    )
                    del level_y[level + 1]

                if level in level_y and level_y[level]:
                    self.renderer.add(
                        self.renderer.v_connection(-8.0 + level, self.t_count,
                                                   level_y[level], 'black')
                    )

            if True:
                if level - 1 in level_y \
                        and level_y[level - 1] == self.t_count - 1:
                    self.renderer.add(
                        self.renderer.v_connection(-8.0 + level, self.t_count,
                                                   self.t_count - 1, 'black')
                    )
                if level + 1 in level_y:
                    del level_y[level + 1]

                if level in level_y and level_y[level]:
                    self.renderer.add(
                        self.renderer.v_connection(-8.0 + level, self.t_count,
                                                   level_y[level], 'black')
                    )


            level_y[level] = self.t_count

            self.render_ensemble_mover_line(ensembles, sub_mp, color='gray')
#            self.render_replica_line(len(ensembles), sub_set, color='gray')

        self.t_count += 1

#        self.render_ensemble_line(ensembles, sset)
#        self.render_replica_line(len(ensembles), sset)

#        self.renderer.add(self.renderer.block(-8.0, self.t_count, 'black'))
#        self.renderer.add(
#            self.renderer.label(-8.0, self.t_count, 3,
#                                'storage.sampleset[%d]' % sset.idx[storage],
#                                align='start', color='black')
#        )

        self.t_count += 1

        self.renderer.shift_x = - 9.0
        self.renderer.shift_y = -0.5
        self.renderer.height = 1.0 * self.t_count + 1.0
        self.renderer.width = 1.0 * len(ensembles) + 20.5

    def render_ensemble_mover_head(self, ensembles, yp=None):
        for ens_idx, ens in enumerate(ensembles):
            txt = chr(ens_idx + 65)

            label = ens.name if hasattr(ens, 'name') else ens.__class__.__name__[:-8]

            self.renderer.add(
                self.renderer.vertical_label(
                    ens_idx, yp, 5, '[' + txt + '] ' + label,
                    align='start'
                )
            )

    def render_ensemble_mover_line(self, ensembles, mover, yp=None, color='black'):
        if yp is None:
            yp = self.t_count

        storage = self.storage
        for ens_idx, ens in enumerate(ensembles):
            txt = chr(ens_idx + 65)
            show = False
            in_ens = mover.input_ensembles
            out_ens = mover.output_ensembles
            if in_ens is None or None in in_ens or ens in in_ens:
                self.renderer.add(
                    self.renderer.connector(ens_idx, yp - 0.12, 'green', ''))

                show = True

            if out_ens is None or None in out_ens or ens in out_ens:
                self.renderer.add(
                    self.renderer.connector(ens_idx, yp + 0.12, 'red', ''))

                show = True

            if show:
                self.renderer.add(
                    self.renderer.block(ens_idx, yp, 'rgb(200,200,200)', txt))


#        for ens_idx, ens in enumerate(ensembles):
#            samp_ens = [samp for samp in sset if samp.ensemble is ens]
#            if len(samp_ens) > 0:
#                traj_idx = samp_ens[0].trajectory.idx[storage]
#                self.ens_x[ens_idx] = traj_idx
#                self.traj_ens_x[traj_idx] = ens_idx
#                self.traj_ens_y[traj_idx] = self.t_count

    def render_ensemble_line(self, ensembles, sset, yp=None, color='black'):
        if yp is None:
            yp = self.t_count

        storage = self.storage
        for ens_idx, ens in enumerate(ensembles):
            samp_ens = [samp for samp in sset if samp.ensemble is ens]
            if len(samp_ens) > 0:
                traj_idx = samp_ens[0].trajectory.idx[storage.trajectories]
                txt = str(traj_idx)
                if len(samp_ens) > 1:
                    txt += '+'

                my_color = color

                if traj_idx in self.ens_x:
                    self.renderer.add(
                        self.renderer.v_hook(self.traj_ens_x[traj_idx],
                                             self.traj_ens_y[traj_idx], ens_idx,
                                             self.t_count, 'black'))
                    if self.traj_ens_x[traj_idx] != ens_idx:
                        my_color = 'red'
                else:
                    if len(self.ens_x) > 0:
                        my_color = 'red'

                if my_color != 'gray':
                    self.renderer.add(
                        self.renderer.connector(ens_idx, yp, my_color, txt))

        for ens_idx, ens in enumerate(ensembles):
            samp_ens = [samp for samp in sset if samp.ensemble is ens]
            if len(samp_ens) > 0:
                traj_idx = samp_ens[0].trajectory.idx[storage.trajectories]
                self.ens_x[ens_idx] = traj_idx
                self.traj_ens_x[traj_idx] = ens_idx
                self.traj_ens_y[traj_idx] = self.t_count

    def render_replica_line(self, replica, sset, yp=None, color='black'):
        if yp is None:
            yp = self.t_count

        storage = self.storage
        for repl_idx in range(0, replica):
            samp_repl = [samp for samp in sset if samp.replica == repl_idx - 1]
            xp = repl_idx + 10
            if len(samp_repl) > 0:
                traj_idx = samp_repl[0].trajectory.idx[storage.trajectories]
                txt = str(traj_idx)
                if len(samp_repl) > 1:
                    txt += '+'
                my_color = color
                if traj_idx in self.repl_x:
                    self.renderer.add(
                        self.renderer.v_hook(self.traj_repl_x[traj_idx],
                                             self.traj_repl_y[traj_idx], xp,
                                             self.t_count, 'black'))
                    if self.traj_repl_x[traj_idx] != xp:
                        my_color = 'red'
                else:
                    if len(self.repl_x) > 0:
                        my_color = 'red'

                if my_color != 'gray':
                    self.renderer.add(
                        self.renderer.connector(xp, yp, my_color, txt))

        for repl_idx in range(replica):
            samp_repl = [samp for samp in sset if samp.replica == repl_idx - 1]
            xp = repl_idx + 10
            if len(samp_repl) > 0:
                traj_idx = samp_repl[0].trajectory.idx[storage.trajectories]
                self.repl_x[repl_idx] = traj_idx
                self.traj_repl_x[traj_idx] = xp
                self.traj_repl_y[traj_idx] = self.t_count

    @staticmethod
    def _get_min_max(d):
        return min(d.values()), max(d.values())


class PathTreeBuilder(object):
    def __init__(self, storage, op=None, states=None):
        self.rejected = False
        self.p_x = dict()
        self.p_y = dict()
        self.obj = list()
        self.storage = storage
        self.renderer = TreeRenderer()
        self.op = op
        self.show_redundant = False
        if states is None:
            states = []
        self.states = states

        self.renderer.style = '''
            .opstree text {
                alignment-baseline: central;
                font-size: 10px;
                text-anchor: middle;
            }
            .opstree text.bw.label {
                text-anchor: end;
            }
            .opstree text.fw.label {
                text-anchor: start;
            }
            .opstree text.block {
                fill: white !important;
                stroke: none !important;
            }
            .opstree .repex {
                fill: blue;
            }
            .opstree .new {
                fill: black;
            }
            .opstree .hop {
                fill: blue;
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
            }
            .opstree line {
                stroke-width: 2px;
            }
            .opstree .label {
                fill: black !important;
            }
            .opstree .connection {
                stroke-dasharray: 3 3;
            }
            .opstree .rejected {
                opacity: 0.3;
            }
            .opstree text {
                font-family: Impact;
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
            '''

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
                    'cls' : ['repex']
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
                    'cls' : ['shooting']
                },
                paths.ForwardShootMover: {
                    'name': 'Shooting',
                    'overlap': 'none',
                    'fw': 'blocks',
                    'bw': 'blocks',
                    'all': 'hidden',
                    'overlap_label': '',
                    'suffix': 'f',
                    'label_position': 'left',
                    'cls' : ['shooting']
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
                    'cls' : ['hop']
                },
                paths.PathReversalMover: {
                    'name': 'hop',
                    'overlap': 'line',
                    'fw': '',
                    'bw': '',
                    'all': 'hidden',
                    'overlap_label': 'Reversal',
                    'suffix': 'h',
                    'label_position': 'left',
                    'cls' : ['reversal']
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
                    'cls' : ['unknown']
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
                    'cls' : ['repex']
                }
            },
            'ui': {
                'trajectory': True,
                'step': True,
                'correlation': True,
                'sample': True
            },
            'settings': {
                'register_rejected': False,
                'time_symmetric': True,
            },
            'geometry': {
                'scale_x' : 12,
                'scale_y' : 24,
                'horizontal_gap' : 0.05
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

    def from_samples(self, samples, clear=True):

        self.renderer.scale_x = self.options['geometry']['scale_x']
        self.renderer.scale_y = self.options['geometry']['scale_y']
        self.renderer.horizontal_gap = self.options['geometry']['horizontal_gap']


        if len(samples) == 0:
            # no samples, nothing to do
            # TODO: Raise an exception or just ignore and don't output anything?
            return

        prev = samples[0].trajectory
        old_tc = 1
        move_list = {}
        step_list = {}
        for step in self.storage.steps:
            for ch in step.change:
                if ch.samples is not None:
                    for trial in ch.samples:
                        step_list[trial] = step
                        move_list[trial] = ch

        p_x = dict()
        p_y = dict()

        if clear:
            self.renderer.clear()

        t_count = 1
        shift = 0

        options = self.options

        x_text_stretch = options['geometry']['scale_y'] / options['geometry']['scale_x']

        assume_reversed_as_same = options['settings']['time_symmetric']

        for sample in samples:
            mover_type = type(sample.mover)
            traj = sample.trajectory

            # get connection to existing
            overlap_reversed = False

            index_bw = None

            for snap_idx in range(len(traj)):
                snap = traj[snap_idx]
                if snap in p_x or assume_reversed_as_same and snap.reversed in p_x:
                    connect_bw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                    index_bw = snap_idx
                    shift_bw = connect_bw - snap_idx
                    break

            new_sample = False
            if index_bw is None:
                index_bw = 0
                index_fw = len(traj) - 1
                # no overlap, so skip
                new_sample = True
                shift = 0
            else:
                for snap_idx in range(len(traj) - 1, -1, -1):
                    snap = traj[snap_idx]
                    if snap in p_x or (assume_reversed_as_same and snap.reversed in p_x):
                        connect_fw = p_x[snap] if snap in p_x else p_x[snap.reversed]
                        index_fw = snap_idx
                        shift_fw = connect_fw - snap_idx
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

            if new_sample:
                view_options = options['mover']['new']
            elif mover_type in options['mover']:
                view_options = options['mover'][mover_type]
            else:
                view_options = options['mover']['unknown']

#            print shift, index_bw, index_fw

            traj_str = str(self.storage.idx(sample.trajectory)) + view_options['suffix']

            cls = [] + view_options['cls']

            if sample in move_list:
                move = step_list[sample].change
                accepted = move.accepted
                if not accepted:
                    cls += ['rejected']

            if view_options['label_position'] == 'left':
                self.renderer.add(
                    self.renderer.label(shift - x_text_stretch, t_count, 1, traj_str, cls=cls + ['left'])
                )
            elif view_options['label_position'] == 'right':
                self.renderer.add(
                    self.renderer.label(shift + len(traj) - 1 + x_text_stretch, t_count, 1, traj_str, cls=cls + ['right'])
                )

            if index_bw > 0:
                self.renderer.add(
                    self.renderer.v_connection(shift + index_bw, p_y[traj[index_bw]], t_count, cls=cls + ['bw', 'connection'])
                )
            if index_fw < len(traj) - 1:
                self.renderer.add(
                    self.renderer.v_connection(shift + index_fw + 1, p_y[traj[index_fw]], t_count, cls=cls + ['fw', 'connection'])
                )

            if view_options['overlap'] == 'line':
                self.renderer.add(
                    self.renderer.range(shift + index_bw, t_count, index_fw - index_bw + 1,
                                        view_options['overlap_label'], cls=cls + ['overlap'] )
                )
                for pos, snapshot in enumerate(sample.trajectory[index_bw:index_fw + 1]):
                    if snapshot not in p_x:
                        p_x[snapshot] = shift + pos + index_bw
                        p_y[snapshot] = t_count

            if view_options['bw'] == 'line':
                self.renderer.add(
                    self.renderer.range(shift + 0, t_count, index_bw, cls=cls + ['bw'])
                )
                for pos, snapshot in enumerate(sample.trajectory[0:index_bw]):
                    p_x[snapshot] = shift + pos
                    p_y[snapshot] = t_count

            if view_options['fw'] == 'line':
                self.renderer.add(
                    self.renderer.range(shift + index_fw + 1, t_count, len(traj) - (index_fw + 1), cls=cls + ['fw'])
                )
                for pos, snapshot in enumerate(sample.trajectory[index_fw+1:]):
                    p_x[snapshot] = shift + pos + index_fw + 1
                    p_y[snapshot] = t_count

            for pos, snapshot in enumerate(sample.trajectory):
                conf_idx = snapshot
                if not conf_idx in p_y or True:
                    pos_x = shift + pos
                    pos_y = t_count

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
                        self.renderer.add(
                            self.renderer.block(
                                pos_x,
                                pos_y,
                                txt,
                                extend_left = pos > 0,
                                extend_right = pos < len(sample.trajectory) - 1,
                                cls = cls + b_cls
                            ))

                        p_x[conf_idx] = pos_x
                        p_y[conf_idx] = pos_y

            t_count += 1

        self.p_x = p_x
        self.p_y = p_y

        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = self._get_min_max(self.p_y)

        self.renderer.shift_x = min_x - 8.5 * x_text_stretch
        self.renderer.shift_y = 2
        self.renderer.height = max_y - min_y + 5.0
        self.renderer.width = max_x - min_x + 10.0 * x_text_stretch

        matrix = self._to_matrix()

        if hasattr(self, 'states') and len(self.states) > 0:
            for color, op in self.states.iteritems():
                xp = None
                for y in range(0, max_y - min_y + 1):
                    left = None
                    yp = y + min_y
                    for x in range(0, (max_x - min_x + 1)):
                        xp = x + min_x

                        # if matrix[y][x] is not None:
                        #     self.renderer.pre(
                        #         self.renderer.shade(xp, yp, 0.9,
                        #                             'black')
                        #     )

                        if matrix[y][x] is not None and bool(op(matrix[y][x])):
                            if left is None:
                                left = xp
                        else:
                            if left is not None:
                                self.renderer.pre(
                                    self.renderer.shade(left, yp, xp - left, cls=[color])
                                )
                                left = None

                    if left is not None:
                        self.renderer.pre(
                            self.renderer.shade(left, yp, xp - left + 1, cls=[color])
                        )

        self.renderer.add(
            self.renderer.label(self.renderer.shift_x + 5.0 * x_text_stretch, 0, 1, 'smp')
        )

        self.renderer.add(
            self.renderer.label(self.renderer.shift_x + 4.0 * x_text_stretch, 0, 1, 'cyc')
        )

        self.renderer.add(
            self.renderer.label(self.renderer.shift_x + 6.5 * x_text_stretch, 0, 1, 'cor')
        )

        for tc, s in enumerate(samples):
            self.renderer.add(
                self.renderer.rect(-10.0 * x_text_stretch + min_x, 1 + tc - 0.45, max_x - min_x + 11 * x_text_stretch, 1, cls=['tableline'])
            )
            if tc > 0 and not paths.Trajectory.is_correlated(s.trajectory, prev):
                self.renderer.add(
                    self.renderer.h_range(self.renderer.shift_x + 6.5 * x_text_stretch, old_tc - 0.1, 1 + tc - old_tc + 0.2, 'black', "" ))

                old_tc = 1 + tc
                prev = s.trajectory

            self.renderer.add(
                self.renderer.label(self.renderer.shift_x + 5.0 * x_text_stretch, 1 + tc, 1, str(
                    self.storage.idx(s)))
            )

            if s in step_list:
                txt = str(step_list[s].mccycle)
            else:
                txt = '---'

            self.renderer.add(
                self.renderer.label(self.renderer.shift_x + 4.0 * x_text_stretch, 1 + tc, 1, str(
                    txt))
            )

        # self.renderer.add(
        #     self.renderer.h_range(self.renderer.shift_x + 2.0, 0.9, len(samples) + 0.2, 'black', "" ))

        self.renderer.add(
            self.renderer.h_range(self.renderer.shift_x + 6.5 * x_text_stretch, old_tc - 0.1, 1 + len(samples) - old_tc + 0.2, 'black', "", extend_bottom=False))

        self.renderer.add(
            self.renderer.label(self.renderer.shift_x + 3.0 * x_text_stretch, 1 + tc, 1, '')
        )

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
     plt.ioff() # turn off interactive mode
     fig=plt.figure(figsize=(2,2))
     ax = fig.add_subplot(111)
     nx.draw_shell(self, ax=ax)
     output = StringIO.StringIO()
     fig.savefig(output,format='svg')
     plt.ion() # turn on interactive mode
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
            nx.draw_networkx_nodes(
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
        data = json_graph.tree_data(self.G,len(self.G)-1)
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
              .text(function(d) { return d.name; });

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
        self.from_samples(self.samples)
        self.view = self.renderer

    def rebuild(self):
        """Rebuild the internal structures.

        It seems like some changes in the visualization require a complete
        rebuild. That's not ideal. If that can be changed, this function
        could be removed.

        """
        self.view.reset_pad()
        self.from_samples(self.samples)


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
                              if s.replica==self.replica]
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
    

import svgwrite
import os

class TreeRenderer(object):
    def __init__(self):
        self.start_x = 0
        self.start_y = 0
        self.scale_x = 24
        self.scale_y = 24
        self.scale_th = 24
        self.font_family = "Futura"
        self.font_size = 0.3
        self.document = None
        self.font_baseline_shift = +0.05
        self.stroke_width = 0.05
        self.horizontal_gap = 0.05
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
        if item in ['block', 'shade', 'v_connection', 'h_connection', 'label']:
            # This will delay the execution of the draw commands until we know where to draw
            return self._delay(object.__getattribute__(self, 'draw_' + item))
        else:
            return object.__getattribute__(self, item)


    def _delay(self, func):
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
        return (self._x(x), self._y(y))

    def _wh(self, w, h):
        return (self._w(w), self._h(h))

    def _xb(self, x, y):
        return (self._x(x), self._y(y + self.font_baseline_shift))


    def _pad(self, xy, wh):
        self.min_x = min(self.min_x, xy[0], xy[0] + wh[0])
        self.min_y = min(self.min_y, xy[1], xy[1] + wh[1])
        self.max_x = max(self.max_x, xy[0] + wh[0], xy[0])
        self.max_y = max(self.max_y, xy[1] + wh[1], xy[1])

    def draw_block(self, x,y, color = "blue", text = "", align="middle", extend_right= True, extend_left = True, padding=None):
        document = self.document
        if padding is None:
            padding = self.horizontal_gap

        self._pad(self._xy(x - 0.5, y - 0.3), self._wh(1.0, 0.6))

        ret = []

        ret.append(document.rect(
            insert = self._xy(x - 0.5 + padding, y - 0.3),
            size = self._wh(1.0 - 2 * padding, 0.6),
            fill = color,
        ))
        if extend_left:
            ret.append(document.circle(
                center = self._xy(x - 0.5, y),
                r = self._w(padding),
                stroke_width = 0,
                stroke = color,
                fill = color
            ))
        if extend_right:
            ret.append(document.circle(
                center = (self._xy(x + 0.5, y)),
                r = self._w( padding),
                stroke_width = 0,
                stroke = color,
                fill = color
            ))
        ret.append(document.text(
            text = str(text)[:4],
            insert = self._xb(x, y),
            text_anchor = align,
            font_size = self._h(self.font_size),
            alignment_baseline = 'middle',
            font_family = self.font_family,
            fill = 'white'
        ))

        return ret

    def draw_shade(self, x, y, w, color, stroke_width = None):
        document = self.document
        if stroke_width is None:
            stroke_width = self.stroke_width

        self._pad(self._xy(x - 0.5, y - 0.35), self._wh(1.0, 0.7))

        return [document.rect(
            insert = self._xy(x - 0.5, y + 0.35),
            size = self._wh(w, 0.1),
            fill = color,
            stroke = color,
            stroke_width = self._th(stroke_width)
        )]

    def draw_shade_alt(self, x, y, w, color, stroke_width = None):
        document = self.document
        if stroke_width is None:
            stroke_width = self.stroke_width

        self._pad(self._xy(x - 0.5, y - 0.35), self._wh(1.0, 0.7))

        return [document.rect(
            insert = self._xy(x - 0.5, y - 0.35),
            size = self._wh(w, 0.7),
            fill = 'none',
            stroke = color,
            stroke_width = self._th(stroke_width)
        )]

    def draw_v_connection(self, x, y1, y2, color, stroke_width = None, padding = None):
        document = self.document
        if stroke_width is None:
            stroke_width = self.stroke_width
        if padding is None:
            padding = self.horizontal_gap

        self._pad(self._xy(x - 0.5 - stroke_width, y1), self._wh(2 * stroke_width, y2-y1))

        return [document.line(
            start = self._xy(x - 0.5, y1 + padding),
            end = self._xy(x - 0.5, y2 - padding),
            stroke_width = self._th(stroke_width),
            stroke = color,
        )]

    def draw_h_connection(self, x1, x2, y, color, stroke_width = None, padding = None):
        document = self.document
        if stroke_width is None:
            stroke_width = self.stroke_width
        if padding is None:
            padding = self.horizontal_gap

        self._pad(self._xy(x1, y - stroke_width), self._wh(x2-x1, 2 * stroke_width))

        return [document.line(
            start = self._xy(x1 + 0.5 + padding, y),
            end = self._xy(x2 - 0.5, y),
            stroke_width = self._th(stroke_width),
            stroke = color,
        )]

    def _text_align_to_int(self, s):
        if s == "start":
            return 1
        elif s == "end":
            return -1
        else:
            return 0

    def draw_label(self, x, y, w, text, align = "middle", color = "black", shift = 0.7):
        document = self.document

        self._pad(
            self._xy(x + shift * self._text_align_to_int(align), y),
            self._wh(w * self._text_align_to_int(align), 0.6)
        )

        return [document.text(
            text = str(text),
            insert = self._xb(x + shift * self._text_align_to_int(align), y),
            text_anchor = align,
            font_size = self._h(self.font_size),
            alignment_baseline = 'middle',
            font_family = self.font_family,
            fill = color
        )]

    def to_document(self):

        self.start_x = self.margin
        self.start_y = self.margin


        self.document = svgwrite.Drawing()
        self.document['width'] = str(self._width()) + 'px'
        self.document['height'] = str(self._height()) + 'px'

        for obj in self.obj:
            parts = obj()
            map(self.document.add, parts)

#        self.document.add(self.document.rect(
#            insert = (0, 0),
#            size = (self._width(), self._height()),
#            fill = 'none',
#            stroke = 'black'
#        ))


        return self.document

    def to_svg(self):
        doc = self.to_document()
        return doc.tostring()

    def to_html(self, svg = None):
        if svg is None:
            svg = self.to_svg()

        html = '<!DOCTYPE html><html style="margin:0px; padding:0px; width:100%;">' + svg + '<body style="margin:0px; padding:0px;"></body></html>'

        return html

    def _height(self):
        return self._h(self.height) + self.margin * 2

    def _width(self):
        return self._w(self.width) + self.margin * 2

    def write_html(self, file = 'tree.html'):
        with open(file, 'w') as f:
            f.write(self.to_html())

    def save_pdf(self, file='tree.pdf'):
        h = self._height()
        w = self._width()

        h_margin = 3.5
        v_margin = 3.5

        page_height =str(25.4 / 75.0 * h + v_margin)+'mm'
        page_width =str(25.4 / 75.0 * w + h_margin)+'mm'

        with open('tree_xxx.html', 'w') as f:
            f.write(self.to_html())

        f.closed

        bashCommand = "wkhtmltopdf -l --page-width " + page_width + " --page-height " + page_height + " --disable-smart-shrinking -B 1mm -L 1mm -R 1mm -T 1mm tree_xxx.html " + file
        os.system(bashCommand)

    def clear(self):
        self.obj = []



class PathTreeBuilder(object):
    def __init__(self, storage, op=None, states = None):
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

    @staticmethod
    def construct_heritage(storage, sample):
        list_of_samples = []

        samp = sample

        while len(samp.details.inputs) > 0:
            if len(samp.details.inputs) == 1:
                # just one sample so use this
                list_of_samples.append(samp)
                samp = samp.details.inputs[0]
            else:
                # if there are more than one input choose the most useful one
                # e.g. for ReplicaExchange the initial one
                found_one = False
                for input in samp.details.inputs:
                    if input.trajectory == list_of_samples[-1].trajectory:
                        # got it
                        found_one = True
                        samp = input
                        break

                if not found_one:
                    break

        # reverse to get origin first
        return [samp for samp in reversed(list_of_samples)]

    def from_samples(self, samples, clear=True):
        if len(samples) == 0:
            # no samples, nothing to do
            # TODO: Raise an exception or just ignore and don't output anything?
            return

        p_x = dict()
        p_y = dict()

        if clear:
            self.renderer.clear()

        t_count = 0
        shift = 0

        lightcolor = "gray"

        for sample in samples:
            if hasattr(sample.details, 'start_point'):
                old_traj = sample.details.start_point.trajectory
                old_index = sample.details.start_point.index
                old_conf = old_traj[old_index].configuration

                new_traj = sample.details.final_point.trajectory
                new_index = sample.details.final_point.index
                new_conf = new_traj[new_index].configuration

                accepted = sample.details.accepted

                if sample.trajectory is new_traj or self.rejected:
                    t_count += 1
                    if not old_conf in p_x:
                        for pos, snapshot in enumerate(old_traj):
                            conf = snapshot.configuration
                            p_x[conf] = pos
                            p_y[conf] = t_count

                            pos_x = p_x[conf]
                            pos_y = p_y[conf]
                            if self.op is not None:
                                self.renderer.add(self.renderer.block(pos_x, pos_y, "black", self.op(snapshot)))
                            else:
                                self.renderer.add(self.renderer.block(pos_x, pos_y, "black", ""))

                        self.renderer.add(
                            self.renderer.label(0, t_count, 1, str(self.storage.idx(new_traj)) + 'b', align='end',color='black')
                        )

                        t_count += 1

                    shift = p_x[old_conf] - new_index


                    fontcolor = "black"

                    draw_okay = False

                    mover_name = ''

                    if hasattr(sample.details, 'mover'):
                        mover_name = sample.details.mover.name

                    if mover_name == "BackwardShootMover":
                        color = "green"
                        if not accepted:
                            color = lightcolor
                            fontcolor = lightcolor
                        self.renderer.add(
                            self.renderer.v_connection(shift + new_index, p_y[old_conf], t_count, color)
                        )
                        self.renderer.add(
                            self.renderer.label(shift, t_count, 1, str(self.storage.idx(new_traj)) + 'b', align='end',color=fontcolor)
                        )
                        draw_okay = True
                    elif mover_name == 'ForwardShootMover':
                        color = "red"
                        if not accepted:
                            color = lightcolor
                            fontcolor = lightcolor
                        self.renderer.add(
                            self.renderer.v_connection(shift + new_index + 1, p_y[old_conf], t_count, color)
                        )
                        self.renderer.add(
                            self.renderer.label(shift + len(new_traj) - 1, t_count, 1, str(self.storage.idx(new_traj)) + 'f', align='start',color=fontcolor)
                        )
                        draw_okay = True

                    if not accepted:
                        color = lightcolor

                    if draw_okay:
                        for pos, snapshot in enumerate(new_traj):
                            conf = snapshot.configuration
                            if not conf in p_y:
                                p_y[conf] = t_count
                                p_x[conf] = shift + pos

                                pos_x = p_x[conf]
                                pos_y = p_y[conf]
                                if self.op is not None:
                                    self.renderer.add(self.renderer.block(pos_x, pos_y, color, self.op(snapshot)))
                                else:
                                    self.renderer.add(self.renderer.block(pos_x, pos_y, color, ""))

        self.p_x = p_x
        self.p_y = p_y

        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = self._get_min_max(self.p_y)

        self.renderer.shift_x = min_x - 1.5
        self.renderer.shift_y = 0
        self.renderer.height = max_y - min_y + 2.0
        self.renderer.width = max_x - min_x + 3.0

        op_names = { arg[0] : arg[1] for arg in self.states }
        ops = {op : self.storage.collectivevariable.load(op) for op in op_names.keys() }

        matrix = self._to_matrix()

        for y in range(0, max_y - min_y + 1):
            rr = { op_name : None for op_name in op_names.keys() }
            yp = y + min_y
            for x in range(0, (max_x - min_x + 1)):
                xp = x + min_x
                for r in rr:
                    op = ops[r]
                    if matrix[y][x] is not None and bool(op(matrix[y][x])):
                        if rr[r] is None:
                            rr[r] = xp
                    else:
                        if rr[r] is not None:
                            self.renderer.pre(
                                self.renderer.shade(rr[r], yp, xp - rr[r], op_names[r])
                            )
                            rr[r] = None

            for r in rr:
                if rr[r] is not None:
                    self.renderer.pre(
                        self.renderer.shade(rr[r], yp, xp - rr[r] + 1, op_names[r])
                    )



    def _get_min_max(self, d):
        return min(d.values()), max(d.values())

    def _to_matrix(self):
        min_x, max_x = self._get_min_max(self.p_x)
        min_y, max_y = self._get_min_max(self.p_y)

        matrix = [[None] * (max_x - min_x + 1) for n in range(max_y - min_y + 1)]

        for s in self.p_x:
            px = self.p_x[s]
            py = self.p_y[s]
            matrix[py - min_y][px - min_x] = s

        return matrix
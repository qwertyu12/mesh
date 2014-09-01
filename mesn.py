
from __future__ import division

import numpy as np
from skimage.util.dtype import dtype_range
from skimage import draw
from skimage import measure

from skimage.viewer.plugins.plotplugin import Plugin


__all__ = ['ConstrainedMeshPlacer']


class ConstrainedMeshPlacer(Plugin):
    """Plugin to compute interpolated intensity under a scan line on an image.

    See PlotPlugin and Plugin classes for additional details.

    Parameters
    ----------
    maxdist : float
        Maximum pixel distance allowed when selecting end point of scan line.
    limits : tuple or {None, 'image', 'dtype'}
        (minimum, maximum) intensity limits for plotted profile. The following
        special values are defined:

            None : rescale based on min/max intensity along selected scan line.
            'image' : fixed scale based on min/max intensity in image.
            'dtype' : fixed scale based on min/max intensity of image dtype.
    """
    name = 'Mesh Placer'

    def __init__(self, maxdist=10, epsilon='deprecated',
                 limits='image', **kwargs):
        super(ConstrainedMeshPlacer, self).__init__(**kwargs)
        self.maxdist = maxdist
#        print(self.help())

    def attach(self, image_viewer):
        super(ConstrainedMeshPlacer, self).attach(image_viewer)
        self.viewer = image_viewer
        image = image_viewer.original_image

        h, w = image.shape[0:2]
        x = [w / 3, 2 * w / 3]
        y = [h / 2] * 2

        num_cells = 10
        xinit = x[0]
        yinit = y[0]        
        
        self.line_tools = [LineToolC(self.image_viewer.ax,maxdist=self.maxdist,on_move=self.line_changed) for xxxxzzzz in range(num_cells*3+1)]
        for i in range(num_cells):
            self.line_tools[i].end_points = np.array([[xinit+(15*i),yinit],[xinit+(15*(i+1)),yinit]])
        for i in range(num_cells,2*num_cells):
            self.line_tools[i].end_points = np.array([[xinit+(15*(i-num_cells)),yinit+15],[xinit+(15*((i-num_cells)+1)),yinit+15]])
        for i in range(num_cells*2,num_cells*3+1):
            step = i - 2*num_cells
            self.line_tools[i].end_points = np.array([[xinit+(15*step),yinit],[xinit+(15*step),yinit+15]])
        
        self.toprange = self.line_tools[:num_cells]
        self.midrange = self.line_tools[num_cells*2:]
        self.bottomrange = self.line_tools[num_cells:num_cells*2]
        
        for i in range(0,len(self.toprange)):
            self.toprange[i].anchors = [self.toprange[0], self.toprange[-1]]
            self.toprange[i].position = 'mid'
            self.bottomrange[i].anchors = [self.bottomrange[0], self.bottomrange[-1]]
            self.bottomrange[i].position = 'mid'
            
        for x in self.midrange[1:-1]:
            x.topanchors = [self.toprange[0], self.toprange[-1]]
            x.bottomanchors = [self.bottomrange[0], self.bottomrange[-1]]
            x.position = 'midspine'
            
        self.toprange[0].position = 'anchorleft'
        self.toprange[-1].position = 'anchorright'
        self.midrange[0].position = 'anchorleftmid'
        self.midrange[-1].position = 'anchorrightmid'
        self.bottomrange[0].position = 'anchorleft'
        self.bottomrange[-1].position = 'anchorright'
        
        for i in [0,-1]:
            self.toprange[i].children = self.toprange
            self.toprange[i].midchildren = self.midrange
            self.bottomrange[i].children = self.bottomrange
            self.bottomrange[i].midchildren = self.midrange
        
        for x in self.line_tools:
            x.ncells = num_cells
        

    def help(self):
        helpstr = ("Line profile tool",
                   "+ and - keys or mouse scroll changes width of scan line.",
                   "Select and drag ends of the scan line to adjust it.")
        return '\n'.join(helpstr)

    def line_changed(self, end_points):
        x, y = np.transpose(end_points)
#        for x in self.line_tools:
#            x.end_points = x.end_points
        
        
        
        
        self.viewer.redraw()
#        for t in self.line_tools:
#            t.redraw()
#        self.line_tool.end_points = end_points

    def output(self):
        """Return the drawn line and the resulting scan.

        Returns
        -------
        line_image : (M, N) uint8 array, same shape as image
            An array of 0s with the scanned line set to 255.
            If the linewidth of the line tool is greater than 1,
            sets the values within the profiled polygon to 128.
        scan : (P,) or (P, 3) array of int or float
            The line scan values across the image.
        """
#        end_points = self.line_tool.end_points
        return [x.end_points for x in self.midrange]

try:
    from matplotlib import lines
except ImportError:
    print("Could not import matplotlib -- skimage.viewer not available.")

from skimage.viewer.canvastools.base import CanvasToolBase, ToolHandles
import math


__all__ = ['LineToolC']


class LineToolC(CanvasToolBase):
    """Widget for line selection in a plot.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Matplotlib axes where tool is displayed.
    on_move : function
        Function called whenever a control handle is moved.
        This function must accept the end points of line as the only argument.
    on_release : function
        Function called whenever the control handle is released.
    on_enter : function
        Function called whenever the "enter" key is pressed.
    maxdist : float
        Maximum pixel distance allowed when selecting control handle.
    line_props : dict
        Properties for :class:`matplotlib.lines.Line2D`.

    Attributes
    ----------
    end_points : 2D array
        End points of line ((x1, y1), (x2, y2)).
    """
    def __init__(self, ax, on_move=None, on_release=None, on_enter=None,
                 maxdist=10, line_props=None):
        super(LineToolC, self).__init__(ax, on_move=on_move, on_enter=on_enter,
                                       on_release=on_release)

        props = dict(color='r', linewidth=1, alpha=0.25, solid_capstyle='butt')
        props.update(line_props if line_props is not None else {})
        self.linewidth = props['linewidth']
        self.maxdist = 10
        self.position = 'blah'
        self._active_pt = None

        x = (0, 0)
        y = (0, 0)
        self._end_pts = np.transpose([x, y])

        self._line = lines.Line2D(x, y, visible=False, animated=True, **props)
        ax.add_line(self._line)

        self._handles = ToolHandles(ax, x, y)
        self._handles.set_visible(False)
        self._artists = [self._line, self._handles.artist]

        if on_enter is None:
            def on_enter(pts):
                x, y = np.transpose(pts)
                print("length = %0.2f" % np.sqrt(np.diff(x)**2 + np.diff(y)**2))
        self.callback_on_enter = on_enter

        self.connect_event('button_press_event', self.on_mouse_press)
        self.connect_event('button_release_event', self.on_mouse_release)
        self.connect_event('motion_notify_event', self.on_move)

    @property
    def end_points(self):
        return self._end_pts

    @end_points.setter
    def end_points(self, pts):
        self._end_pts = np.asarray(pts)
        self._line.set_data(np.transpose(pts))
        self._handles.set_data(np.transpose(pts))
        self._line.set_linewidth(self.linewidth)

        self.set_visible(True)
        self.redraw()

    def on_mouse_press(self, event):
        if event.button != 1 or not self.ax.in_axes(event):
            return
        self.set_visible(True)
        idx, px_dist = self._handles.closest(event.x, event.y)
        if px_dist < self.maxdist:
            self._active_pt = idx
#        else:
#            self._active_pt = 0
#            x, y = event.xdata, event.ydata
#            self._end_pts = np.array([[x, y], [x, y]])
    
    def distalong(self, anchors, d):
        x1,y1 = anchors[0].end_points[0]
        x2,y2 = anchors[1].end_points[1]
        m = (y2-y1)/(x2-x1)
        xadj = x1 + math.sqrt((d*d)/(m*m+1))
        yadj = m*(xadj - x1)+y1
        
        return [xadj, yadj]

    def distance(self, xy1, xy2):
        x1,y1 = xy1
        x2,y2 = xy2

        return math.sqrt((y1-y2)**2 + (x1-x2)**2)        
    
    def closest(self, p1, anchors):
        x3,y3 = p1
        x1,y1 = anchors[0].end_points[0]
        x2,y2 = anchors[1].end_points[1]
        
#        print x1,y1,'\n',x2,y2,'\n',x3,y3        
        
        m = (y2-y1)/(x2-x1)
        xadj = (x3+m*(y3-y1)+m*m*x1)/(m*m+1)
        yadj = m*(xadj - x1)+y1
        return [xadj, yadj]

    def interpolate(self, points):
        if self.position == 'anchorright':
            return [self.closest(points[0], self.anchors),points[1]]
        elif self.position == 'anchorleft':
            return [points[0],self.closest(points[1], self.anchors)]
        elif self.position == 'mid':
            return [self.closest(x, self.anchors) for x in points]
        elif self.position == 'midspine':
            return [self.closest(points[0], self.topanchors), self.closest(points[1], self.bottomanchors)]

        return points

    def on_mouse_release(self, event):
        if event.button != 1:
            return
        self._active_pt = None
        self.callback_on_release(self.geometry)

    def on_move(self, event):
#        print 'onmove'
        idx, px_dist = self._handles.closest(event.x, event.y)
#        if px_dist < self.maxdist:
#            print px_dist, self.maxdist
        if event.button != 1 or self._active_pt is None:
            return
        if not self.ax.in_axes(event):
            return
        self.update(event.xdata, event.ydata)
        self.callback_on_move(self.geometry)

    def update(self, x=None, y=None, orig=None):
        if x is not None:
            self._end_pts[self._active_pt, :] = x, y
        
        if orig is None:
            px = np.round(self._end_pts)
            px = self.interpolate(px)
            self._end_pts = np.asarray(px)            
            self.end_points = self._end_pts
        
        
        if orig == None and self.position == 'anchorleft' or self.position == 'anchorright':
            fulld = self.distance(self.anchors[0].end_points[0],self.anchors[1].end_points[1])
            step = fulld / self.ncells
            
            topa = self.midchildren[1].topanchors
            bota = self.midchildren[1].bottomanchors
            
            topd = self.distance(topa[0].end_points[0], topa[1].end_points[1])
            botd = self.distance(bota[0].end_points[0], bota[1].end_points[1])
            
            topstep = topd/self.ncells
            botstep = botd/self.ncells

            for i,x in enumerate(self.midchildren[1:-1]):
                x.end_points = [x.distalong(x.topanchors,topstep*(i+1)), x.distalong(x.bottomanchors, botstep*(i+1))]
            
            if self.position == 'anchorleft' and self._active_pt == 0:
                self.end_points = [self.end_points[0], self.distalong(self.anchors,step)]
                for i,x in enumerate(self.children):
                    if x is not self:
                        if x.position != 'anchorright':
                            x.end_points = [x.distalong(x.anchors,step*i), x.distalong(x.anchors, step*(i+1))]
                        else:
                            x.end_points = [x.distalong(x.anchors,step*i), x.end_points[1]]
#                        x.update(orig = i)
            if self.position == 'anchorright' and self._active_pt == 1:
                self.end_points = [self.distalong(self.anchors,step), self.end_points[1]]
                for i,x in enumerate(self.children):
                    if x is not self:
                        if x.position != 'anchorleft':
                            x.end_points = [x.distalong(x.anchors,step*i), x.distalong(x.anchors, step*(i+1))]
                        else:
                            x.end_points = [x.end_points[0], x.distalong(x.anchors,step)]
#                        x.update(orig = i)

    @property
    def geometry(self):
        return self.end_points

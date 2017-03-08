#!/usr/bin/env python2
# encoding: utf-8

"""\
Judge forward-folded candidates in computational protein design pipelines.

Usage:
    show_my_designs.py [options] <pdb_directories>...
    show_my_designs.py --version

Options:
    -F, --no-fork
        Do not fork into a background process.

    -f, --force
        Force the cache to be regenerated.

    -q, --quiet
        Build the cache, but don't launch the GUI.

    -v, --version
        Print the version number and exit.

Features:
    1. Extract quality metrics from forward-folded models and plot them against 
       each other in any combination.

    2. Easily visualize specific models by right-clicking on plotted points.  
       Add your own visualizations by writing `*.sho' scripts.

    3. Plot multiple designs at once, for comparison purposes.

    4. Keep notes on each design, and search your notes to find the designs you 
       want to visualize.
"""

## Imports
import collections, glob, gzip, os, re, shutil, subprocess, sys
import gtk, pango, yaml
import matplotlib.pyplot as plt, numpy as np, scipy as sp, pandas as pd

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg


class Design (object):

    def __init__(self, directory, use_cache=True):
        self.directory = directory
        self.cache_path = os.path.join(directory, 'models.pkl')
        self.notes_path = os.path.join(directory, 'notes.txt')
        self.rep_path = os.path.join(directory, 'representative.txt')

        self._models = None
        self._notes = ""
        self._representative = None

        self._load_models(use_cache)
        self._load_annotations()

    def __str__(self):
        return '<ModelGroup dir={}>'.format(self.directory)

    def __len__(self):
        return len(self.paths)


    @property
    def paths(self):
        return self._models['path']

    @property
    def notes(self):
        return self._notes

    @notes.setter
    def notes(self, notes):
        self._notes = notes
        self._save_notes()

    @property
    def representative(self):
        if self._representative is None:
            return np.argmin(self.get_metric('total_score'))
        else:
            return self._representative

    @representative.setter
    def representative(self, index):
        self._representative = index
        self._save_representative()

    @property
    def representative_path(self):
        return self.paths[self.representative]

    def get_metric(self, metric):
        if metric not in self.defined_metrics:
            message = "No such metric: '{}'\n".format(metric)
            message += "Defined metrics are: " + ', '.join(
                    "'{}'".format(x) for x in self.defined_metrics)
            print type(metric), ' '.join(str(type(x)) for x in self.defined_metrics)

            raise RuntimeError(message)

        return self._models[metric]

    def get_coord(self, x_metric, y_metric, index=None):
        i = index if index is not None else self.representative
        return self.get_metric(x_metric)[i], self.get_metric(y_metric)[i]

    @property
    def defined_metrics(self):
        # This method returns the names of any column in self._models that 
        # contains numeric data.  Any dtype other than 'object' is assumed to 
        # be numeric.
        return {x for x in self._models.keys()
                if self._models[x].dtype != 'object'}


    def _load_models(self, use_cache):
        """
        Load a variety of score and distance metrics for the structures found 
        in the given directory.  As much information as possible will be 
        cached.  Note that new information will only be calculated for file 
        names that haven't been seen before.  If a file changes or is deleted, 
        the cache will not be updated to reflect this and you may be presented 
        with stale data.
        """

        # Make sure the given directory matches all of our expectations: i.e.  
        # that it exists and contains PDB files.

        if not os.path.exists(self.directory):
            raise IOError("'{}' does not exist".format(self.directory))
        if not os.path.isdir(self.directory):
            raise IOError("'{}' is not a directory".format(self.directory))
        if not os.listdir(self.directory):
            raise IOError("'{}' is empty".format(self.directory))
        if not glob.glob(os.path.join(self.directory, '*.pdb*')):
            raise IOError("'{}' doesn't contain any PDB files".format(self.directory))

        # Find all the structures in the given directory, then decide which 
        # have already been cached and which haven't.

        pdb_paths = glob.glob(os.path.join(self.directory, '*.pdb*'))
        base_pdb_names = set(os.path.basename(x) for x in pdb_paths)

        if use_cache and os.path.exists(self.cache_path):
            cached_records = pd.read_pickle(self.cache_path).to_dict('records')
            cached_paths = set(
                    record['path'] for record in cached_records
                    if 'path' in record)
            uncached_paths = [
                    pdb_path for pdb_path in pdb_paths
                    if os.path.basename(pdb_path) not in cached_paths]
        else:
            cached_records = []
            uncached_paths = pdb_paths

        # Calculate score and distance metrics for the uncached paths, then 
        # combine the cached and uncached data into a single data frame.

        uncached_records = parse_records_from_pdbs(uncached_paths)
        self._models = pd.DataFrame(cached_records + uncached_records)

        # Make sure at least two metrics have been associated with each model 
        # in this directory.

        if len(self.defined_metrics) == 0:
            raise IOError("no metrics defined for the models in '{}'".format(self.directory))
        if len(self.defined_metrics) == 1:
            defined_metric = self.defined_metrics.pop()
            raise IOError("only found one metric '{}' for the models in '{}', need at least two".format(defined_metric, self.directory))

        # If everything else looks good, cache the data frame so we can load 
        # faster next time.

        if not self._models.empty:
            self._models.to_pickle(self.cache_path)

    def _load_annotations(self):
        try:
            with open(self.notes_path) as file:
                self._notes = file.read()
        except IOError:
            pass

        try:
            with open(self.rep_path) as file:
                self._representative = int(file.read())
        except IOError:
            pass

    def _save_notes(self):
        with open(self.notes_path, 'w') as file:
            file.write(self.notes)

        if os.path.exists(self.notes_path) and not self.notes:
            os.remove(self.notes_path)

    def _save_representative(self):
        if self._representative is not None:
            with open(self.rep_path, 'w') as file:
                file.write(str(self._representative))

        elif os.path.exists(self.rep_path):
            os.remove(self.rep_path)


class ShowMyDesigns (gtk.Window):

    def __init__(self, designs):
        
        # Setup the parent class.

        gtk.Window.__init__(self)
        self.add_events(gtk.gdk.KEY_PRESS_MASK)
        self.connect('key-press-event', self.on_hotkey_press)
        self.set_icon_from_file(os.path.join(
            os.path.dirname(os.path.realpath(__file__)), 'icon.png'))

        # Setup the data members.

        self.designs = designs
        self.keys = list()
        self.selected_model = None

        self.defined_metrics = sorted(set.intersection(
                *[x.defined_metrics for x in self]))
        self.x_metric = (
                default_x_metric
                if default_x_metric in self.defined_metrics
                else self.defined_metrics[1])
        self.y_metric = (
                default_y_metric
                if default_y_metric in self.defined_metrics
                else self.defined_metrics[1])

        # Setup the GUI.

        self.connect('destroy', lambda x: gtk.main_quit())
        self.set_default_size(int(1.618 * 529), 529)

        model_list = self.setup_model_list()
        model_viewer = self.setup_model_viewer()
        #model_viewer.set_size_request(529, 529)

        if len(designs) > 1:
            widgets = gtk.HPaned()
            widgets.add1(model_list)
            widgets.add2(model_viewer)
        else:
            widgets = model_viewer

        self.add(widgets)
        self.set_border_width(3)
        self.update_everything()
        self.show_all()
        self.set_focus(None)

    def __iter__(self):
        return iter(self.designs.values())


    def setup_model_list(self):
        list_store = gtk.ListStore(str)

        text = gtk.CellRendererText()
        icon = gtk.CellRendererPixbuf()

        self.view = gtk.TreeView(list_store)
        self.view.set_model(list_store)
        self.view.set_rubber_banding(True)
        self.view.set_enable_search(False)
        self.view.set_headers_visible(False)

        columns = [
                ('Name', 'directory'),
        ]

        for index, parameters in enumerate(columns):
            title, attr = parameters

            def cell_data_func(column, cell, model, iter, attr):
                key = model.get_value(iter, 0)
                design = self.designs[key]
                text = getattr(design, attr)
                cell.set_property('text', text)

            def sort_func(model, iter_1, iter_2, attr):
                key_1 = model.get_value(iter_1, 0)
                key_2 = model.get_value(iter_2, 0)
                design_1 = self.designs[key_1]
                design_2 = self.designs[key_2]
                value_1 = getattr(design_1, attr)
                value_2 = getattr(design_2, attr)
                return cmp(value_1, value_2)

            list_store.set_sort_func(index, sort_func, attr);

            column = gtk.TreeViewColumn(title, text)
            column.set_cell_data_func(text, cell_data_func, attr)
            column.set_sort_column_id(index)
            self.view.append_column(column)

        selector = self.view.get_selection()
        selector.connect("changed", self.on_select_designs)
        selector.set_mode(gtk.SELECTION_MULTIPLE)

        scroller = gtk.ScrolledWindow()
        scroller.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        scroller.add(self.view)

        frame = gtk.Frame()
        frame.add(scroller)

        self.search_form = gtk.Entry()
        self.search_form.set_icon_from_stock(gtk.ENTRY_ICON_SECONDARY, gtk.STOCK_FIND)

        search_buffer = self.search_form.get_buffer()
        search_buffer.connect('deleted-text', self.on_search_in_notes)
        search_buffer.connect('inserted-text', self.on_search_in_notes)

        save_paths = gtk.Button("Save selected paths")
        save_funnels = gtk.Button("Save selected funnels")

        save_paths.connect(
                'clicked', lambda _: self.save_interesting_paths())
        save_funnels.connect(
                'clicked', lambda _: self.save_interesting_funnels())

        vbox = gtk.VBox()
        vbox.pack_start(self.search_form, expand=False)
        vbox.pack_start(frame)
        vbox.pack_start(save_paths, expand=False)
        vbox.pack_start(save_funnels, expand=False)

        return vbox

    def setup_model_viewer(self):
        plot = self.setup_plot()
        notes = self.setup_annotation_area()

        panes = gtk.VPaned()
        panes.add1(plot)
        panes.add2(notes)

        return panes

    def setup_plot(self):
        figure = Figure(facecolor='#edecea')

        # Create the axes.

        self.axes = figure.add_axes((0.15, 0.15, 0.75, 0.75))
        self.axes.set_ylabel('Score')

        # Create the canvas.

        self.canvas = FigureCanvas(figure)
        self.canvas.mpl_connect('pick_event', self.on_select_model)
        self.canvas.mpl_connect('button_press_event', self.on_click_plot_mpl)
        self.canvas.mpl_connect('motion_notify_event', self.on_move_mouse_mpl)
        self.canvas.connect('button-press-event', self.on_click_plot_gtk)
        self.canvas.set_size_request(-1, 350)

        # Create the tool bar.

        self.toolbar = NavigationToolbar(self.canvas, self)

        # Create the axis menus.
        
        x_axis_menu = gtk.combo_box_new_text()
        y_axis_menu = gtk.combo_box_new_text()

        for i, metric in enumerate(self.defined_metrics):
            x_axis_menu.append_text(get_metric_title(metric))
            y_axis_menu.append_text(get_metric_title(metric))

            if metric == self.x_metric: x_axis_menu.set_active(i)
            if metric == self.y_metric: y_axis_menu.set_active(i)

        x_axis_menu.connect('changed', self.on_change_x_metric)
        y_axis_menu.connect('changed', self.on_change_y_metric)

        # Place all the widgets.

        self.mouse_position = gtk.Label("")

        table = gtk.Table(3, 5)
        table.attach(self.toolbar, 0, 1, 0, 3)
        table.attach(self.mouse_position, 3, 4, 1, 2, xoptions=0, yoptions=0, xpadding=3)

        vbox = gtk.VBox()
        vbox.pack_start(self.canvas)
        vbox.pack_start(table, expand=False)

        return vbox

    def setup_annotation_area(self):
        self.notes = gtk.TextView()
        self.notes.set_wrap_mode(gtk.WRAP_WORD)
        self.notes.set_size_request(-1, 100)
        self.notes.set_left_margin(3)
        self.notes.set_right_margin(3)
        self.notes.set_pixels_above_lines(3)
        self.notes.set_pixels_below_lines(3)
        self.notes.set_cursor_visible(True)
        self.notes.get_buffer().connect('changed', self.on_edit_annotation)

        scroll_window = gtk.ScrolledWindow()
        scroll_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        scroll_window.add(self.notes)

        frame = gtk.Frame()
        frame.add(scroll_window)

        return frame


    def on_hotkey_press(self, widget, event):
        key = gtk.gdk.keyval_name(event.keyval).lower()
        if event.state & gtk.gdk.CONTROL_MASK: key = 'ctrl-' + key
        if event.state & gtk.gdk.SHIFT_MASK: key = 'shift-' + key
    
        hotkeys = {
                'escape': self.normal_mode,
        }
        normal_mode_hotkeys = {
                'i': self.insert_mode,     'a': self.insert_mode,
                'z': self.zoom_mode,
                'x': self.pan_mode,
                'c': self.refocus_plot,
                'tab': self.cycle_y_metric,
                'space': self.cycle_x_metric,
        }
        multi_design_hotkeys = {
                'j': self.next_design,     'f': self.next_design,
                'k': self.previous_design, 'd': self.previous_design,
                'slash': self.search_mode,
        }

        if self.get_focus() not in (self.notes, self.search_form):
            hotkeys.update(normal_mode_hotkeys)
            if len(self.designs) > 1:
                hotkeys.update(multi_design_hotkeys)

        if key in hotkeys:
            hotkeys[key]()
            return True

    def on_search_in_notes(self, entry_buffer, *_):
        self.update_filter()

    def on_select_designs(self, selection) :
        new_keys = []
        old_keys = self.keys[:]
        self.keys = []
        model, paths = selection.get_selected_rows()

        for path in paths:
            iter = model.get_iter(path)
            key = model.get_value(iter, 0)
            new_keys.append(key)

        # Don't change the order of designs that were already selected.  The 
        # order affects how the color of the design in the score vs rmsd plot, 
        # and things get confusing if it changes.

        for key in old_keys:
            if key in new_keys:
                self.keys.append(key)

        for key in new_keys:
            if key not in self.keys:
                self.keys.append(key)

        # Rename the window based on the current selection.

        subtitle = ""
        if len(self.keys) == 1:
            subtitle = " ({})".format(self.keys[0])
        if len(self.keys) > 1:
            subtitle = " ({}, ...)".format(self.keys[0])

        self.set_title("Show My Designs" + subtitle)

        # This is an efficiency thing.  The 'J' and 'K' hotkeys works in two 
        # steps: first unselect everything and then select the next row in 
        # order.  Redrawing the plot is expensive, so it's worthwhile to skip 
        # redrawing after that first step.

        if self.keys:
            self.update_plot()
            self.update_annotations()

    def on_select_model(self, event):
        self.selected_model = event.ind[0], event.artist.design

    def on_move_mouse_mpl(self, event):
        if event.xdata is None or event.ydata is None:
            # The data coordinates will be None only if the mouse is outside 
            # the data area.
            self.mouse_position.set_text("")
        else:
            coord = '{:0.2f}, {:0.2f}'.format(event.xdata, event.ydata)
            self.mouse_position.set_text(coord)

    def on_click_plot_mpl(self, event):
        pass

    def on_click_plot_gtk(self, widget, event):
        # Ignore any event that isn't a right button click or a left button 
        # click with the control key held down.

        is_right_click = \
                (event.button == 3) or \
                (event.button == 1 and event.get_state() & gtk.gdk.CONTROL_MASK)

        if not is_right_click: return
        if self.toolbar._active == 'PAN': return
        if self.toolbar._active == 'ZOOM': return
        if self.selected_model is None: return

        # Figure out which model was clicked.

        index, design = self.selected_model
        path = os.path.join(design.directory, design.paths[index])
        is_rep = (design.representative == index)
        self.selected_model = None

        # Search for scripts that can perform some action using the clicked 
        # model.  Such scripts must have the `*.sho' suffix and may be located 
        # anywhere from the directory containing the models to any directory 
        # below that.  Any scripts that are found will be used to populate a 
        # drop-down menu.  If selected, the script will be called with sh as 
        # the interpreter and the path to the model as the singular argument.

        directory = os.path.abspath(design.directory)
        sho_scripts = []

        while directory != os.path.abspath('/'):
            sho_pattern = os.path.join(directory, '*.sho')
            sho_scripts += glob.glob(sho_pattern)
            directory = os.path.dirname(directory)

        # Create and display the drop-down menu.

        file_menu = gtk.Menu()

        for script in sho_scripts:
            title = os.path.basename(os.path.splitext(script)[0])
            title = title[0].upper() + title[1:]
            title = title.replace('_', ' ')

            item = gtk.MenuItem(title)
            item.connect('activate',
                    lambda *args: try_to_run_command([script, path]))
            file_menu.append(item)

        view_in_pymol = gtk.MenuItem("View model in pymol")
        view_in_pymol.connect('activate',
                lambda *args: try_to_run_command(['pymol', path]))
        file_menu.append(view_in_pymol)

        view_in_chimera = gtk.MenuItem("View model in chimera")
        view_in_chimera.connect('activate',
                lambda *args: try_to_run_command(['chimera', path]))
        file_menu.append(view_in_chimera)

        file_menu.append(gtk.SeparatorMenuItem())

        copy_path = gtk.MenuItem("Copy path to model")
        copy_path.connect('activate', self.on_copy_model_path, path)
        file_menu.append(copy_path)

        if index == design.representative:
            choose_rep = gtk.MenuItem("Reset representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, None)
        else:
            choose_rep = gtk.MenuItem("Set as representative")
            choose_rep.connect(
                'activate', self.on_set_representative, design, index)
        file_menu.append(choose_rep)

        file_menu.foreach(lambda item: item.show())
        file_menu.popup(None, None, None, event.button, event.time)

    def on_copy_model_path(self, widget, path):
        import subprocess
        xsel = subprocess.Popen(['xsel', '-pi'], stdin=subprocess.PIPE)
        xsel.communicate(path)

    def on_set_representative(self, widget, design, index):
        design.set_representative(index)
        self.update_plot()

    def on_edit_annotation(self, buffer):
        assert len(self.keys) == 1
        design = self.designs[self.keys[0]]
        bounds = buffer.get_bounds()
        design.notes = buffer.get_text(*bounds)

    def on_change_x_metric(self, widget):
        self.x_metric = widget.get_active_text()
        self.update_plot()

    def on_change_y_metric(self, widget):
        self.y_metric = widget.get_active_text()
        self.update_plot()


    def normal_mode(self):
        self.set_focus(None)

        if self.toolbar._active == 'PAN':
            self.toolbar.pan()

        if self.toolbar._active == 'ZOOM':
            self.toolbar.zoom()

    def insert_mode(self):
        self.set_focus(self.notes)

    def search_mode(self):
        self.set_focus(self.search_form)

    def zoom_mode(self):
        self.toolbar.zoom()

    def pan_mode(self):
        self.toolbar.pan()

    def refocus_plot(self):
        self.toolbar.home()
        self.normal_mode()

    def next_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        num_paths = model.iter_n_children(None)
        if paths[-1][0] < model.iter_n_children(None) - 1:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[-1][0] + 1)
            self.view.scroll_to_cell(paths[-1][0] + 1)

    def previous_design(self):
        selection = self.view.get_selection()
        model, paths = selection.get_selected_rows()
        if paths[0][0] > 0:
            for path in paths: selection.unselect_path(path)
            selection.select_path(paths[0][0] - 1)
            self.view.scroll_to_cell(paths[0][0] - 1)

    def cycle_x_metric(self):
        i = self.defined_metrics.index(self.x_metric)
        i = (i + 1) % len(self.defined_metrics)
        if self.defined_metrics[i] == self.y_metric:
            i = (i + 1) % len(self.defined_metrics)
        
        # Change the axis by programmatically selecting a new entry in the 
        # corresponding drop-down menu in the toolbar.  This is incredibly 
        # roundabout (and it kinda breaks encapsulation, although I consider 
        # ModelViewer and ModelToolbar to be friends), but it's the only way I 
        # know to keep the drop-down menu in sync.

        self.toolbar.x_axis_menu.set_active(i)

    def cycle_y_metric(self):
        i = self.defined_metrics.index(self.y_metric)
        i = (i + 1) % len(self.defined_metrics)
        if self.defined_metrics[i] == self.x_metric:
            i = (i + 1) % len(self.defined_metrics)

        # Change the axis by programmatically selecting a new entry in the 
        # corresponding drop-down menu in the toolbar.  This is incredibly 
        # roundabout (and it kinda breaks encapsulation, although I consider 
        # ModelViewer and ModelToolbar to be friends), but it's the only way I 
        # know to keep the drop-down menu in sync.

        self.toolbar.y_axis_menu.set_active(i)

    def save_interesting_paths(self):
        chooser = gtk.FileChooserDialog(
                "Save interesting paths",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('interesting_paths.txt')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            selected_designs = [self.designs[key] for key in self.keys]
            with open(chooser.get_filename(), 'w') as file:
                file.writelines(
                        os.path.join(
                            design.directory,
                            design.paths[design.representative]) + '\n'
                        for design in selected_designs)

        chooser.destroy()

    def save_interesting_funnels(self):
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt

        selected_designs = [self.designs[key] for key in self.keys]

        chooser = gtk.FileChooserDialog(
                "Save interesting funnels",
                parent=self,
                action=gtk.FILE_CHOOSER_ACTION_SAVE,
                buttons=(
                    gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                    gtk.STOCK_SAVE, gtk.RESPONSE_OK))

        chooser.set_current_folder(os.getcwd())
        chooser.set_current_name('interesting_funnels.pdf')

        response = chooser.run()

        if response == gtk.RESPONSE_OK:
            pdf = PdfPages(chooser.get_filename())

            for index, design in enumerate(selected_designs):
                plt.figure(figsize=(8.5, 11))
                plt.suptitle(design.directory)

                self.plot_models(plt.gca(), [design])

                pdf.savefig()
                plt.close()

            pdf.close()

        chooser.destroy()

    def plot_models(self, axes, designs, **kwargs):
        from itertools import count

        labels = kwargs.get('labels', None)
        x_metric = kwargs.get('x_metric', self.x_metric)
        y_metric = kwargs.get('y_metric', self.y_metric)

        # Define the colors that the plot will use.

        red =    '#ef2929', '#cc0000', '#a40000'
        orange = '#fcaf3e', '#f57900', '#ce5c00'
        yellow = '#fce94f', '#edd400', '#c4a000'
        green =  '#8ae234', '#73d216', '#4e9a06'
        blue =   '#729fcf', '#3465a4', '#204a87'
        purple = '#ad7fa8', '#75507b', '#5c3566'
        brown =  '#e9b96e', '#c17d11', '#8f5902'

        def color_from_cycle(index):
            cycle = (blue[1], red[1], green[2], orange[1], purple[1], brown[1],
                     blue[0], red[0], green[1], orange[0], purple[0], brown[0])
            return cycle[index % len(cycle)]

        # Clear the axes and reset the axis labels

        axes.clear()
        axes.set_xlabel(get_metric_title(x_metric))
        axes.set_ylabel(get_metric_title(y_metric))

        # Plot the two axes.

        for index, design in enumerate(designs):
            rep = design.representative
            x = design.get_metric(x_metric)
            y = design.get_metric(y_metric)
            color = color_from_cycle(index)
            label = labels[index] if labels is not None else ''
            size = np.clip(7500 / (len(x)), 2, 15)

            # Highlight the representative model.

            axes.scatter(
                    [x[rep]], [y[rep]],
                    s=60, c=yellow[1], marker='o', edgecolor='none')

            # Draw the whole score vs distance plot.

            lines = axes.scatter(
                    x, y,
                    s=size, c=color, marker='o', edgecolor='none',
                    label=label, picker=True)

            lines.paths = design.paths
            lines.design = design

        # Pick the axis limits based on the range of every design.  This is done 
        # so you can scroll though every design without the axes changing size.

        def get_metric_limits(metric):
            values = np.concatenate([x.get_metric(metric) for x in self])
            try: return metric_limits[metric](values)
            except KeyError: return min(values), max(values)

        x_min, x_max = get_metric_limits(x_metric)
        y_min, y_max = get_metric_limits(y_metric)

        x_pad = 0.05 * (x_max - x_min)
        y_pad = 0.05 * (y_max - y_min)

        axes.set_ylim(
            bottom=y_min - y_pad,
            top=y_max + y_pad,
        )
        axes.set_xlim(
            left=x_min - x_pad,
            right=x_max + x_pad,
        )

        # If appropriate, draw guides for the given axes.

        x_guide = metric_guides.get(self.x_metric)
        y_guide = metric_guides.get(self.y_metric)

        if x_guide is not None:
            axes.axvline(x_guide, color='gray', linestyle='--')
        if y_guide is not None:
            axes.axhline(y_guide, color='gray', linestyle='--')

        # If between 1 and 5 designs are being shown, show a legend.

        if labels and 1 < len(designs) < 5:
            axes.legend(loc='upper right')


    def update_everything(self):
        self.update_filter()
        self.update_annotations()
        self.update_plot()

    def update_plot(self):
        designs = [self.designs[k] for k in self.keys]
        self.plot_models(self.axes, designs, labels=self.keys)
        self.canvas.draw()

    def update_annotations(self):
        if len(self.keys) == 1:
            design = self.designs[self.keys[0]]
            self.notes.get_buffer().set_text(design.notes)
            self.notes.set_sensitive(True)
        else:
            self.notes.set_sensitive(False)
        

    def update_filter(self):
        model = self.view.get_model()
        selector = self.view.get_selection()
        model.clear()

        def query_matches_design(design):
            needle = self.search_form.get_text()
            haystack = design.notes

            if needle.islower():
                haystack = haystack.lower()

            return needle in haystack

        for key in sorted(self.designs):
            if query_matches_design(self.designs[key]):
                model.append([key])

        selector.select_path((0,))


class FigureCanvas (FigureCanvasGTKAgg):

    def __init__(self, figure):
        FigureCanvasGTKAgg.__init__(self, figure)

    def button_press_event(self, widget, event):
        FigureCanvasGTKAgg.button_press_event(self, widget, event)
        return False


class NavigationToolbar (NavigationToolbar2GTKAgg):

    toolitems = ( # (fold)
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        #(None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        #(None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
    )

    def __init__(self, canvas, parent):
        NavigationToolbar2GTKAgg.__init__(self, canvas, parent)

        def make_axis_menu(initial_metric, callback):
            store = gtk.ListStore(str, str)
            combo_box = gtk.ComboBox(store)
            cell = gtk.CellRendererText()
            combo_box.pack_start(cell, True)
            combo_box.add_attribute(cell, 'text', 1)

            for i, metric in enumerate(parent.defined_metrics):
                store.append([metric, get_metric_title(metric)])
                if metric == initial_metric:
                    combo_box.set_active(i)

            combo_box.connect('changed', callback)
            return combo_box

        self.x_axis_menu = make_axis_menu(
                parent.x_metric, parent.on_change_x_metric)
        self.y_axis_menu = make_axis_menu(
                parent.y_metric, parent.on_change_y_metric)

        table = gtk.Table(3, 4)
        table.attach(gtk.SeparatorToolItem(), 0, 1, 0, 3)
        table.attach(self.y_axis_menu, 1, 2, 1, 2, xoptions=0, yoptions=0)
        table.attach(gtk.Label(' vs. '), 2, 3, 1, 2, xoptions=0, yoptions=0)
        table.attach(self.x_axis_menu, 3, 4, 1, 2, xoptions=0, yoptions=0)

        tool_item = gtk.ToolItem()
        tool_item.add(table)

        self.insert(tool_item, len(self.toolitems))

    def set_message(self, message):
        pass



metric_titles = {
        'total_score': 'Total Score (REU)',
        'loop_rmsd': u'Loop RMSD (Å)',
        'delta_buried_unsats': u'Δ Buried Unsats',
}

metric_limits = {
        'total_score': lambda x: (
            min(x),
            np.percentile(x, 85)),

        'loop_rmsd': lambda x: (
            0.025 * max(x),
            max(x)),
}
metric_guides = {
        'loop_rmsd': 1.0,
}


default_x_metric = 'loop_rmsd'
default_y_metric = 'total_score'

def show_my_designs(directories, use_cache=True, launch_gui=True, fork_gui=True):
    try:
        try:
            designs = load_designs(directories, use_cache=use_cache)
        except IOError as error:
            if str(error):
                print "Error:", str(error)
                sys.exit()
            else:
                raise

        if designs and launch_gui:
            # If the user wants to run in a background process, try to fork.  
            # But for some reason fork() doesn't seem to work on Macs, so just 
            # run the GUI in the main process if anything goes wrong.
            try:
                if fork_gui and os.fork():
                    sys.exit()
            except Exception:
                pass

            gui = ShowMyDesigns(designs)
            gtk.main()

    except KeyboardInterrupt:
        print

def load_designs(directories, use_cache=True):
    designs = collections.OrderedDict()

    for directory in directories:
        designs[directory] = Design(directory, use_cache)

    return designs

def parse_records_from_pdbs(pdb_paths):
    records = []

    for i, path in enumerate(pdb_paths):

        # Update the user on our progress, because this is often slow.

        sys.stdout.write("\rReading '{}' [{}/{}]".format(
            os.path.dirname(path), i+1, len(pdb_paths)))
        sys.stdout.flush()

        # Read the PDB file, which we are assuming is gzipped.

        try:
            def smart_open(path):
                if path.endswith('.gz'): return gzip.open(path)
                else: return open(path)

            with smart_open(path) as file:
                lines = file.readlines()

        except IOError:
            print "\nFailed to read '{}'".format(path)
            continue

        # Parse the pdb file.  This method may be overloaded to parse 
        # different kinds of information.

        record = {'path': os.path.basename(path)}
        parse_record_from_pdb(record, path, lines)
        records.append(record)

    if pdb_paths: print
    return records

def parse_record_from_pdb(record, pdb_path, lines):
    # Get different information from different lines in the PDB file.  Some 
    # of these lines are specific to certain simulations.

    for line in lines:
        if line.startswith('total_score') or line.startswith('pose'):
            record['total_score'] = float(line.split()[1])

        if line.startswith('loop_backbone_rmsd'):
            record['loop_rmsd'] = float(line.split()[1])

        if line.startswith('delta_buried_unsats'):
            record['delta_buried_unsats'] = float(line.split()[1])

def try_to_run_command(command):
    with open(os.devnull, 'w') as devnull:
        try: subprocess.Popen(command, stdout=devnull)
        except OSError as error:
            message = gtk.MessageDialog(
                    parent=None, 
                    flags=0, 
                    type=gtk.MESSAGE_ERROR, 
                    buttons=gtk.BUTTONS_OK,
            )
            message.set_markup("<b>Failed to run {}</b>".format(command[0]))
            message.format_secondary_text(str(error))
            message.run()
            message.destroy()

def get_metric_title(metric):
    naive_title = metric.title()
    naive_title = naive_title.replace('_', ' ')
    naive_title = naive_title.replace('-', ' ')
    return metric_titles.get(metric, naive_title)


def main():
    import docopt
    args = docopt.docopt(__doc__)

    if args['--version']:
        from show_my_designs import __version__
        print '{0} {1} (python {2[0]}.{2[1]})'.format(
                os.path.basename(sys.argv[0]), __version__, sys.version_info)
        raise SystemExit

    show_my_designs(
            args['<pdb_directories>'],
            use_cache=not args['--force'],
            launch_gui=not args['--quiet'],
            fork_gui=not args['--no-fork'],
    )


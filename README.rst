Show My Designs
===============
The purpose of this program is to make it easier to judge forward-folded 
candidates in computational protein design pipelines.

Most computational design pipelines have one step that searches for sequences 
which satisfy the design goals and another that simulates some of those 
sequences to see which really do satisfy those goals.  This second step is 
usually called forward-folding or computational validation.  For each design 
chosen to be validated, hundreds of forward-folding simulations may be run. 
Each of these produces a single model which can be characterized by metrics 
expressing how realistic it is and how well it satisfies the design goals. 
These metrics include force-field scores, RMSDs to target structures, buried 
unsatisfied H-bond counts, and possibly other things like that.  A design is 
promising if the forward-folded models that best satisfy the design goals are 
also the most realistic.

This program provides a number of utilities and features to make it easier to 
find promising designs:

1. Extract quality metrics from forward-folded models and plot them against 
   each other in any combination.

2. Easily visualize specific models by right-clicking on plotted points.

3. Plot multiple designs at once, for comparison purposes.

4. Keep notes on each design, and search your notes to find the designs you 
   want to visualize.

Installation
------------
The most difficult part of installing show_my_designs is making sure all the 
required dependencies are present:

- python 2.7
- pygtk
- numpy
- scipy
- pandas
- numexpr

On linux, these most of these should already be installed, and all of these 
should be available through whatever package manager your distribution uses. 
On mac, homebrew seems to be the most promising route, but I haven't tried it.

Bugs and New Features
---------------------
If you find a bug, open an issue through the github interface::

    https://github.com/Kortemme-Lab/show_my_designs/issues

If you'd like to fix a bug or make an improvement to the code, fork the project 
and make a pull request::

    https://github.com/Kortemme-Lab/show_my_designs/fork

Usage
-----
Use the `-h` flag to get help on using ``show_my_designs``:

    $ ./show_my_designs.py -h
    Usage: ...

Generally, the only arguments you need are the names of one or more directories 
containing the forward-folded models in the PDB format.  For example::

    $ ls
    show_my_designs.py
    design_1/
    design_2/
    ...

    $ ls design_1
    model_1.pdb
    model_2.pdb
    ...

    $ ./show_my_designs.py design_*

This last command will launch the GUI.  If you specified more than one design 
on the command line, the GUI will have a panel on the left listing all the 
designs being compared.  You can control what is plotted by selecting one or 
more designs from this list.  The search bar at the top of this panel can be 
used to filter the list for designs that have the search term in their 
descriptions.  The buttons at the bottom can be used to save information about 
whatever designs are selected.  The "Save selected paths" button will save a 
text file listing the path to the lowest scoring model for each selected 
design.  The "Save selected funnels" button will save a PDF with the plot for 
each selected design on a separate page.

The upper right area of the GUI will contain a plot with different metrics on 
the two axes where each point represents a single model.  You can right-click 
on any point to take an action on the model represented by that point.  Usually 
this means visualizing the model in an external program, like pymol or chimera. 
You can also run your own custom scripts; see the "customization" section below 
for more information.  

The tool bar below the plot can be used to pan around, zoom in or out, save an 
image of the plot, or change the axes.  If the mouse is over the plot, its 
coordinates will be shown just to the right of these controls.  Below the plot 
is a text form which can be used to enter a description of the design.  These 
descriptions can be searched.  I like using the '+', '++', ... convention to 
rank designs so I can easily search for increasingly good designs.

Customization
-------------
Because every protein design pipeline is different, ``show_my_designs`` was 
written to be flexible.  Providing a new way to visualize specific models is 
trivial. You just need to write a script with the extension ``*.sho`` that 
takes the path of a model as its only argument.  ``show_my_designs`` will 
search for scripts with this extension in every directory starting with the 
directory containing the model in question and going down all the way to the 
root of the file system. Any scripts that are found are added to the menu you 
get by right-clicking on a point, using simple rules (the first letter is 
capitalized and underscores are converted to spaces) to convert the file name 
into a menu item name.

Another common modification is to change what metrics are extracted for each 
model.  By default, only the metrics that outputted by rosetta's loop modeling 
framework are extracted.  The metrics include: the rosetta fullatom score, the 
RMSD to the native backbone, and the number of buried unsatisfied H-bonds.  To 
add new metrics, you have to monkey-patch the ``show_my_designs`` module and 
call ``show_my_designs.main()`` from your own script.

This is more clear with an example.  Say your forward-folding simulation 
outputs an auxiliary file for each model containing all sorts of metrics 
relevant to your particular system.  You can add support for these metrics by 
reimplementing ``show_my_designs.parse_records_from_pdbs()``. This function 
takes a list of paths to PDB files that haven't been cached yet and returns a 
list containing ``{'metric_name': metric_value}`` dictionaries for each one.  
The information in this list is cached so that it doesn't have to be 
regenerated unless necessary.

If your custom metrics are encoded in the PDB file itself, you can reimplement 
``show_my_designs.parse_record_from_pdb()`` instead.  This function is called 
by ``show_my_designs.parse_records_from_pdbs()`` and with a list of the lines 
in a specific PDB file.  It is expected to return the ``{'metric_name': 
metric_value}`` dictionary for that model.

Hotkeys
-------
- j,f,down:   Select the next design, if there is one.
- k,d,up:     Select the previous design, if there is one.
- i,a:        Focus on the description form.
- z:          Use the mouse to zoom on a rectangle.
- x:          Use the mouse to pan (left-click) or zoom (right-click).
- c:          Return to the original plot view.
- slash:      Focus on the search bar.
- tab:        Change the y-axis metric.
- space:      Change the x-axis metric.
- escape:     Unfocus the search and description forms.

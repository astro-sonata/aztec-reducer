# sonata-reducer

Reduction pipelines for SONATA. Currently the only one developed
is for Bok.

### Installation
To install for running the pipeline run the following commands
```
git clone https://github.com/noahfranz13/sonatapy.git
cd sonatapy
pip install .
```

### Running the Pipeline
When you pip install the package a terminal command is automatically
installed as well. The only pipeline that currently exists is for
Bok and can be run as follows
```
reduce-bok -d /path/to/data/directory/ -s standard_star_name -n number_of_cores
```

To run a more detailed reduction of the data you can import to a
jupyter notebook using
```
from sonatapy import Bok
```

### Developer Instructions

To install for development
```
git clone https://github.com/noahfranz13/sonatapy.git
cd sonatapy
pip install -e .
```
So that your pip installation is editable. The code is formatted as
follows:

1. All of the pip installable code is in `py/sonatpy`
2. `telescope.py` is a superclass for any telescope that SONATA uses
   and has the primary reduction methods.
3. `bok.py` is a subclass of `telescope.py` and adds some properties
   of the Bok telescope as instance variables and a reduction pipeline
   function.
4. `exceptions.py` are a few custom exceptions used throughout the code
5. `_version.py` is the one place that holds the code version. This is
   really only important once we have pypi releases.
6. `pipeline` holds the pipeline scripts. If you add a new pipeline for
   a different telescope of different type of data from Bok you have to
   open the `pyproject.toml` and under the `project.scripts` heading add
   `cmd = path.to.command:main` where you replace `cmd` with whatever
   you want the command name to be. An example for Bok is shown.

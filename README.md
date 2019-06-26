# HTInspect

# Description

Tools for inspecting the results of HiTIME searches of liquid-chromatography mass spectrometry data.

# License

HTInspect is released under the terms of the 3-Clause BSD License. The terms of the License are provided in the file called `LIC
ENSE` in the top level of the repository.

# Requirements

HTInspect is written in Python 3.7 and requires the Qt5 toolkit. 

# Installation

1. Clone this repository:

```
git clone https://github.com/mgleeming/HTInspect
```

2. Create a Python virtual environment on your computer, and install HTInspect into the sandbox:
```
cd HTInspect
python3 -m venv sandbox
source sandbox/bin/activate
pip install -U .
```

# Usage

With the sandbox from above activated, run the `HTInspect` program:

```
HTInspect
```

If you do not have the sandbox activated, it can be reactivated like so:
```
cd HTInspect
source sandbox/bin/activate
```

# Temper.py

Try it live at: https://sintel.pythonanywhere.com/

## What is this?

A library for implementing various kinds of routines that are used in the study of [regular temperaments](https://en.wikipedia.org/wiki/Regular_temperament), a niche topic in microtonal music theory. This mostly comes down to 'integer linear algebra', such as Hermite normal forms, solving linear diophantine equations, lattice problems, etc.

It also implements a small web app using flask, where (expert) users can query and search properties for their systems of interest.

## Using the library

See [example.py](example.py) for an example on how to do your own calculations using the library.
Currently it is not distributed on pypi so you will just have to clone it.

## Running the web app locally

The easiest way to run is by just making a venv and installing Flask in there.
You need at least Python 3.10.
Check out the tutorial here: https://flask.palletsprojects.com/en/1.1.x/installation/.
After activating the venv, run `pip install -r requirements.txt`.
Then start the app using `python -u flask_app.py`.

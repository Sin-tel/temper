import traceback
import flask as f
import os

# workaround for git import errors
os.environ["GIT_PYTHON_REFRESH"] = "quiet"
os.environ["GIT_PYTHON_GIT_EXECUTABLE"] = "/usr/bin/git"
import git
from timeout import time_limit, TimeoutException
from info import *
from search import temperament_search

# from werkzeug.middleware.profiler import ProfilerMiddleware


app = f.Flask(__name__)
# app.wsgi_app = ProfilerMiddleware(app.wsgi_app)


@app.route("/index")
@app.route("/")
def index():
    if f.request.method == "GET":
        args = f.request.args
        return f.render_template("./temper.jinja", args=args)
    return f.render_template("./temper.jinja")


@app.route("/search", methods=["GET"])
def search():
    if f.request.method == "GET":
        args = f.request.args.to_dict()
        if len(args) == 0:
            return f.render_template("./search.jinja")

        for k in list(args.keys()):
            if args[k] == "" or args[k].isspace():
                del args[k]

        try:
            with time_limit(5):
                res = temperament_search(args)
        except TimeoutException as e:
            raise TimeoutException("Calculation took too long!") from e

        if f.request.headers.get("X-Requested-With") == "XMLHttpRequest":
            return f.render_template("./search_results.jinja", res=res)

        # Otherwise return the full page
        return f.render_template("./search.jinja", args=args, res=res)

    raise ValueError("nothing submitted")


@app.route("/result", methods=["GET"])
def result():
    if f.request.method == "GET":
        args = f.request.args

        args = args.to_dict()
        for k in list(args.keys()):
            if args[k] == "" or args[k].isspace():
                del args[k]

        args["tenney"] = "tenney" in args

        if "subgroup" in args:
            basis, s_expanded = parse_subgroup(args["subgroup"])
        else:
            # infer basis from commas
            assert (
                "submit_comma" in args and "commas" in args
            ), "subgroup can only be inferred from commas"
            s_expanded = p_limit(97)
            basis = np.eye(len(s_expanded), dtype=np.int64)
            commas = parse_intervals(args["commas"], basis, s_expanded)
            commas = np.hstack(commas)
            subset = [not np.all(k == 0) for k in commas]
            s_expanded = [prime for (prime, s) in zip(s_expanded, subset) if s]
            basis = np.eye(len(s_expanded), dtype=np.int64)

        if "submit_edo" in args:
            assert "edos" in args, "no edos entered"
            T, T_expanded = from_edos(args["edos"], basis, s_expanded)
        elif "submit_comma" in args:
            assert "commas" in args, "no commas entered"
            T, T_expanded = from_commas(args["commas"], basis, s_expanded)
        else:
            raise ValueError("nothing submitted")

        try:
            with time_limit(5):
                html_info = info(T, basis, T_expanded, s_expanded, args)
        except TimeoutException as e:
            raise TimeoutException("Calculation took too long!") from e

        return f.render_template("./result.jinja", res=html_info, args=f.request.args)
    raise ValueError("nothing submitted")


@app.route("/update", methods=["POST"])
def update():
    if f.request.method == "POST":
        wd = os.getcwd()
        if wd.endswith("temper"):
            repo = git.Repo(".")
        else:
            repo = git.Repo("./temper")

        origin = repo.remotes.origin
        if not "main" in repo.heads:
            repo.create_head("main", origin.refs.main)
        repo.heads.main.set_tracking_branch(origin.refs.main).checkout()
        origin.pull()
        return "", 200
    return "", 400


@app.route("/xen")
def xen():
    return '<img src="static/xen.png" width="550">'


@app.route("/test")
def test():
    return "hey5"


@app.errorhandler(500)
def internal_error(_):
    print("500 error caught")
    return "<pre>" + traceback.format_exc() + "</pre>", 500


if __name__ == "__main__":
    # app.run(host="127.0.0.1", debug=True)
    app.run(debug=True)

    # app.run(debug=True, threaded=True)
    # app.run(threaded=True)

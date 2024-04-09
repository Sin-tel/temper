import flask as f
import traceback
import git
from timeout import time_limit, TimeoutException
from info import *
from markupsafe import Markup
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
            if args[k] == "":
                del args[k]

        res = temperament_search(args)
        return f.render_template("./search.jinja", args=args, res=res)

    raise ValueError("nothing submitted")


@app.route("/result", methods=["GET"])
def result():
    if f.request.method == "GET":
        args = f.request.args

        args = args.to_dict()
        for k in list(args.keys()):
            if args[k] == "":
                del args[k]

        args["tenney"] = "tenney" in args

        basis, s_expanded = parse_subgroup(args["subgroup"])
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
    return "hey"


@app.errorhandler(500)
def internal_error(exception):
    print("500 error caught")
    return "<pre>" + traceback.format_exc() + "</pre>"


if __name__ == "__main__":
    # app.run(host="127.0.0.1", debug=True)
    app.run(debug=True)

    # app.run(debug=True, threaded=True)
    # app.run(threaded=True)

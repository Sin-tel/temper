import flask as f
from util import *
import traceback
import git
from timeout import time_limit, TimeoutException


app = f.Flask(__name__)


@app.route("/")
def index():
    return f.render_template("./temper.html")


@app.route("/result", methods=["POST", "GET"])
def result():
    if f.request.method == "GET":
        args = f.request.args

        args = args.to_dict()

        args["tenney"] = "tenney" in args
        args["reduce"] = args.get("reduce")

        if "submit_edo" in args:
            temp = from_edos(args)
        elif "submit_comma" in args:
            temp = from_commas(args)
        else:
            return "error"

        try:
            with time_limit(5):
                html_info = info(temp, args)
        except TimeoutException as e:
            raise TimeoutException("Calculation took too long!") from e

        print("", flush=True)
        return f.render_template("./result.html", res=html_info)


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
    else:
        return "", 400


@app.route("/test")
def test():
    return "succes!\n v0.1.14"


@app.errorhandler(500)
def internal_error(exception):
    print("500 error caught")
    return "<pre>" + traceback.format_exc() + "</pre>"


if __name__ == "__main__":
    app.run(debug=True)
    # app.run()

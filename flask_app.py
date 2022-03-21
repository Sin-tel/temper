import flask as f
from util import *
import git

app = f.Flask(__name__)


@app.route('/')
def index():
	return f.render_template("./temper.html")


@app.route('/result', methods=['POST', 'GET'])
def result():
	if f.request.method == 'GET':
		args = f.request.args

		options = dict()

		options["tenney"] = ("tenney" in args)
		options["reduce"] = ("reduce" in args)

		if "submit_edo" in args:
			temp = from_edos(args)
		elif "submit_comma" in args:
			temp = from_commas(args)
		else:
			return "error"

		html_info = info(temp, options)

		print("", flush=True)
		return f.render_template("./result.html", res=html_info)


@app.route('/update', methods=['POST'])
def update():
	if f.request.method == 'POST':
		repo = git.Repo('./temper')
		origin = repo.remotes.origin
		if not 'main' in repo.heads:
			repo.create_head('main', origin.refs.main)
		repo.heads.main.set_tracking_branch(origin.refs.main).checkout()
		origin.pull()
		return '', 200
	else:
		return '', 400


@app.route('/test')
def test():
	return "succes!\n v0.1.9"


if __name__ == '__main__':
	app.run(debug=True)

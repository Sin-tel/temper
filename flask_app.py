import flask as f
from main import *
import git

app = f.Flask(__name__)

@app.route('/')
def index():
   return f.render_template("./temper.html")

@app.route('/result',methods = ['POST', 'GET'])
def result():
   if f.request.method == 'GET':
      args = f.request.args

      temp = from_edos(args)
      html_info = info(temp)
      return f.render_template("./result.html",res = html_info)

@app.route('/update', methods=['POST'])
def update():
   if f.request.method == 'POST':
      repo = git.Repo('./temper')
      origin = repo.remotes.origin
      repo.create_head('main', origin.refs.main).set_tracking_branch(origin.refs.main).checkout()
      origin.pull()
      return '', 200
   else:
      return '', 400

@app.route('/test')
def test():
   return "succes5!"


if __name__ == '__main__':
   app.run(debug = True)


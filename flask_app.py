import flask as f
from main import *
import subprocess

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

@app.route('/pull')
def pull():
   return subprocess.check_output("/home/sintel/temper/pull.sh")

@app.route('/test')
def test():
   return "succes!"


if __name__ == '__main__':
   app.run(debug = True)


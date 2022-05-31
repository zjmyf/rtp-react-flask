from flask import Flask, render_template
from markupsafe import escape
from flask import url_for
from requests import request

app = Flask(__name__)


@app.route('/hello/')
@app.route('/hello/<name>')
def hello(name=None):
    return render_template('hello.html', name=name)


@app.route('/login')
def login():
    return 'login'


@app.route('/user/<username>')
def profile(username):
    return f'{username}\'s profile'


with app.test_request_context():
    # print(url_for('index'))
    print(url_for('login'))
    print(url_for('login', next='/'))
    print(url_for('profile', username='John Doe'))


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/predict/<smi>')
def predict(smi):
    #long_prb, middle_prob
    # pred_res=ml_pred(smi)
    return {'req': smi}

#searchword = request.args.get('key', '')


if __name__ == '__main__':
    app.run()

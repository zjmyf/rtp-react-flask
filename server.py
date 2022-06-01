from flask import Flask, render_template
from markupsafe import escape
from flask import url_for
from flask import request
from ml.smi_pred import smi_pred

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/predict', methods=['GET', 'POST'])
def predict():
    smi = request.json.get('smi', '')
    return smi_pred(smi)

#searchword = request.args.get('key', '')


if __name__ == '__main__':
    app.run()

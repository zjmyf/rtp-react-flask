from flask import Flask, render_template
from markupsafe import escape
from flask import url_for
from flask import request

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/predict')
def predict():
    #long_prb, middle_prob
    # pred_res=ml_pred(smi)
    smi = request.args.get('smi', '')
    return {'req': smi}

#searchword = request.args.get('key', '')


if __name__ == '__main__':
    app.run()

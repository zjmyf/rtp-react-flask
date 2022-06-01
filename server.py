from flask import Flask, render_template
from flask import request
from flask import Response
from flask import send_file
from ml.smi_pred import smi_pred
import os

from rdkit_script.smi2img import smi2img

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/predict', methods=['GET', 'POST'])
def predict():
    smi = request.json.get('smi', '')
    return smi_pred(smi)


@app.route('/smiVis', methods=['GET', 'POST'])
def smiVis():
    smi = request.json.get('smi', '')
    randId = request.json.get('randId', '')
    print(randId)
    fileStr = smi2img(smi,randId)
    return send_file('static/media/smiImg/'+fileStr,mimetype='image/png')
    
#searchword = request.args.get('key', '')

@app.route('/cleanImg', methods=['GET', 'POST'])
def cleanImg():
    randId = request.json.get('randId', '')
    os.remove('static/media/smiImg/{}.png'.format(randId))
    return {'clean':'ok'}

if __name__ == '__main__':
    app.run()

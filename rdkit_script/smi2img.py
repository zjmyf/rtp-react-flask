from genericpath import exists
from importlib.resources import path
from random import random
from rdkit import Chem
from rdkit.Chem import Draw
from datetime import datetime
import os


def smi2img(smi,randId):
    mol = Chem.MolFromSmiles(smi)
    if (not os.path.exists('static/media/smiImg')):
        os.mkdir('static/media/smiImg')
    
    Draw.MolToFile(mol, 'static/media/smiImg/{}.png'.format(randId))
    return '{}.png'.format(randId)

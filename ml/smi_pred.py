import pickle
from ml.gen_des import smi2fp
from rdkit import Chem

with open('model/best_model.pkl', 'rb') as f:
    best_clf = pickle.load(f)

with open('model/best_scaler.pkl', 'rb') as f:
    std_scaler = pickle.load(f)


def smi_pred(smi):
    dic = {'res': '', 'prob': 0}
    screen_smi_long_set = set()
    screen_smi_middle_set = set()

    mol = Chem.MolFromSmiles(smi)
    if (mol is not None):
        fp_res = smi2fp(smi, 'morgan+')
        fp_res = std_scaler.transform(fp_res)

        pred = best_clf.predict(fp_res)
        long_prob, middle_prob, s_prob = best_clf.predict_proba(fp_res)[0]
        if (pred[0] == 'long'):
            dic['res'] = 'long'
            dic['prob'] = long_prob
        elif (pred[0] == 'middle'):
            dic['res'] = 'middle'
            dic['prob'] = middle_prob
        else:
            dic['res'] = 'too short/none'
            dic['prob'] = s_prob

    else:
        dic['res'] = 'invalid SMILES!'
        dic['prob'] = -1

    return dic

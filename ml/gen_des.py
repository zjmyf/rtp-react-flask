from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import numpy as np
import pandas as pd


class Quantity_des:
    des_dic = {}
    tag_method_dic = {
        'BalabanJ': Descriptors.BalabanJ,
        'BertzCT': Descriptors.BertzCT,
        #'Ipc': Descriptors.Ipc,
        'HallKierAlpha': Descriptors.HallKierAlpha,
        'Kappa1': Descriptors.Kappa1,
        'Kappa2': Descriptors.Kappa2,
        'Kappa3': Descriptors.Kappa3,
        'Chi0': Descriptors.Chi0,
        'Chi1': Descriptors.Chi1,
        'Chi0n': Descriptors.Chi0n,
        'Chi1n': Descriptors.Chi1n,
        'Chi2n': Descriptors.Chi2n,
        'Chi3n': Descriptors.Chi3n,
        'Chi4n': Descriptors.Chi4n,
        'Chi0v': Descriptors.Chi0v,
        'Chi1v': Descriptors.Chi1v,
        'Chi2v': Descriptors.Chi2v,
        'Chi3v': Descriptors.Chi3v,
        'Chi4v': Descriptors.Chi4v,
        'MolLogP': Descriptors.MolLogP,
        'MolMR': Descriptors.MolMR,
        'MolWt': Descriptors.MolWt,
        'ExactMolWt': Descriptors.ExactMolWt,
        'HeavyAtomCount': Descriptors.HeavyAtomCount,
        'HeavyAtomMolWt': Descriptors.HeavyAtomMolWt,
        'NHOHCount': Descriptors.NHOHCount,
        'NOCount': Descriptors.NOCount,
        'NumHAcceptors': Descriptors.NumHAcceptors,
        'NumHDonors': Descriptors.NumHDonors,
        'NumHeteroatoms': Descriptors.NumHeteroatoms,
        'NumRotatableBonds': Descriptors.NumRotatableBonds,
        'NumValenceElectrons': Descriptors.NumValenceElectrons,
        'NumAromaticRings': Descriptors.NumAromaticRings,
        'NumSaturatedRings': Descriptors.NumSaturatedRings,
        'NumAliphaticRings': Descriptors.NumAliphaticRings,
        'FractionCSP3': Descriptors.FractionCSP3,
        'NumSpiroAtoms': Chem.rdMolDescriptors.CalcNumSpiroAtoms,
        'NumBridgeheadAtoms': Chem.rdMolDescriptors.CalcNumBridgeheadAtoms,
        'TPSA': Descriptors.TPSA,
        'LabuteASA': Descriptors.LabuteASA,
        'PEOE_VSA1': Descriptors.PEOE_VSA1,
        'PEOE_VSA2': Descriptors.PEOE_VSA2,
        'PEOE_VSA3': Descriptors.PEOE_VSA3,
        'PEOE_VSA4': Descriptors.PEOE_VSA4,
        'PEOE_VSA5': Descriptors.PEOE_VSA5,
        'PEOE_VSA6': Descriptors.PEOE_VSA6,
        'PEOE_VSA7': Descriptors.PEOE_VSA7,
        'PEOE_VSA8': Descriptors.PEOE_VSA8,
        'PEOE_VSA9': Descriptors.PEOE_VSA9,
        'PEOE_VSA10': Descriptors.PEOE_VSA10,
        'PEOE_VSA11': Descriptors.PEOE_VSA11,
        'PEOE_VSA12': Descriptors.PEOE_VSA12,
        'PEOE_VSA13': Descriptors.PEOE_VSA13,
        'PEOE_VSA14': Descriptors.PEOE_VSA14,
        'SMR_VSA1': Descriptors.SMR_VSA1,
        'SMR_VSA2': Descriptors.SMR_VSA2,
        'SMR_VSA3': Descriptors.SMR_VSA3,
        'SMR_VSA4': Descriptors.SMR_VSA4,
        'SMR_VSA5': Descriptors.SMR_VSA5,
        'SMR_VSA6': Descriptors.SMR_VSA6,
        'SMR_VSA7': Descriptors.SMR_VSA7,
        'SMR_VSA8': Descriptors.SMR_VSA8,
        'SMR_VSA9': Descriptors.SMR_VSA9,
        'SMR_VSA10': Descriptors.SMR_VSA10,
        'SlogP_VSA1': Descriptors.SlogP_VSA1,
        'SlogP_VSA2': Descriptors.SlogP_VSA2,
        'SlogP_VSA3': Descriptors.SlogP_VSA3,
        'SlogP_VSA4': Descriptors.SlogP_VSA4,
        'SlogP_VSA5': Descriptors.SlogP_VSA5,
        'SlogP_VSA6': Descriptors.SlogP_VSA6,
        'SlogP_VSA7': Descriptors.SlogP_VSA7,
        'SlogP_VSA8': Descriptors.SlogP_VSA8,
        'SlogP_VSA9': Descriptors.SlogP_VSA9,
        'SlogP_VSA10': Descriptors.SlogP_VSA10,
        'SlogP_VSA11': Descriptors.SlogP_VSA11,
        'SlogP_VSA12': Descriptors.SlogP_VSA12,
        'EState_VSA1': Descriptors.EState_VSA1,
        'EState_VSA2': Descriptors.EState_VSA2,
        'EState_VSA3': Descriptors.EState_VSA3,
        'EState_VSA4': Descriptors.EState_VSA4,
        'EState_VSA5': Descriptors.EState_VSA5,
        'EState_VSA6': Descriptors.EState_VSA6,
        'EState_VSA7': Descriptors.EState_VSA7,
        'EState_VSA8': Descriptors.EState_VSA8,
        'EState_VSA9': Descriptors.EState_VSA9,
        'EState_VSA10': Descriptors.EState_VSA10,
        'EState_VSA11': Descriptors.EState_VSA11,
    }

    def __init__(self, mol_df):
        self.mol_df = mol_df
        self.smi_list = list(mol_df['smiles'])

    def compute_all(self):
        for smi in self.smi_list:
            self.gen_des(smi)

    def get_res(self):
        self.des_dic={}
        self.compute_all()
        res_df = pd.DataFrame(self.des_dic)
        res_df['class'] = self.mol_df['class']
        res_df['smiles'] = self.mol_df['smiles']
        return res_df

    def gen_des(self, smi):
        mol = Chem.MolFromSmiles(smi)
        for tag, method in self.tag_method_dic.items():
            temp = method(mol)
            self.des_dic.setdefault(tag, []).append(temp)

class Finger_des:
    fp_list=[]
    
    def __init__(self, mol_df):
        self.mol_df = mol_df
        self.smi_list = list(mol_df['smiles'])
    
    def gen_finger(self,smi,finger):
        mol = Chem.MolFromSmiles(smi)
        if (finger=='morgan'):
            fp=AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            self.fp_list.append(list(fp))
            
        elif (finger=='maccs'):
            fp = AllChem.GetMACCSKeysFingerprint(mol)
            self.fp_list.append(list(fp)[1:])
            
        elif (finger=='atom_pairs'):
            fp = AllChem.GetHashedAtomPairFingerprintAsBitVect(mol)
            self.fp_list.append(list(fp))
            
        elif (finger=='rdkit'):
            fp = Chem.RDKFingerprint(mol)
            self.fp_list.append(list(fp))
    
    def get_res(self,finger):
        self.fp_list=[]
        for smi in self.smi_list:
            self.gen_finger(smi,finger)
            
        fp_array = np.array(self.fp_list)
        if (finger == 'maccs'):
            fp_df = pd.DataFrame(
                fp_array, columns=[finger + '-' + str(n) for n in range(166)])
        else:
            fp_df = pd.DataFrame(
                fp_array, columns=[finger + '-' + str(n) for n in range(2048)])

        fp_df['class'] = self.mol_df['class']
        fp_df['smiles'] = self.mol_df['smiles']
        
        return fp_df

    def compute_all(self):
        for smi in self.smi_list:
            self.gen_des(smi)


class Combined_des:
    combined_list = []

    def __init__(self, mol_df):
        self.mol_df = mol_df
        self.smi_list = list(mol_df['smiles'])

    def get_res(self, finger, quanti=False):     
        finger_des = Finger_des(self.mol_df)
        finger_df = finger_des.get_res(finger)
        
        maccs_des = Finger_des(self.mol_df)
        maccs_df = maccs_des.get_res('maccs')
        
        combined_df = pd.merge(finger_df,maccs_df,how='inner',on=['smiles','class'])
        if (quanti==False):
            return combined_df
        else:
            quan_des = Quantity_des(self.mol_df)
            quan_df = quan_des.get_res()
            combined_df2 = pd.merge(combined_df,quan_df,how='inner',on=['smiles','class'])
            return combined_df2

def getDes(mol_df,finger):
    if (finger in ['morgan','maccs','atom_pairs','rdkit']):
        finger_des=Finger_des(mol_df)
        return finger_des.get_res(finger)
    elif (finger=='quanti'):
        quanti=Quantity_des(mol_df)
        return quanti.get_res()
    elif (finger in ['morgan+','atom_pairs+','rdkit+']):
        comb=Combined_des(mol_df)
        return comb.get_res(finger[:-1])
    elif (finger in ['morgan++','atom_pairs++','rdkit++']):
        comb=Combined_des(mol_df)
        return comb.get_res(finger[:-2],quanti=True)
    

def smi2fp(smi,finger):
    df=pd.DataFrame([[smi,'too short/none']],columns=['smiles','class'])
    res=getDes(df,finger)
    return res.drop(['smiles','class'],axis=1)
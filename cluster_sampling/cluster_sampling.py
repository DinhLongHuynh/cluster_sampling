from sklearn.metrics import silhouette_score
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import cluster
from rdkit import Chem
from tqdm import tqdm
from rdkit.Chem import Descriptors
from sklearn.preprocessing import StandardScaler
import warnings

class Cluster_Sampling:
    def __init__(self, df, smiles_column):
        self.smiles_column = smiles_column
        self.df = df[[smiles_column]]
        other_columns = [col for col in df.columns if col != smiles_column]
        if other_columns:
            warnings.warn("Columns other than the SMILES column will be dropped: " + ", ".join(other_columns))
            

        

    def standardize_smiles(self, smiles,desalt,remove_stereo):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        try:
            if desalt:
                frags = Chem.GetMolFrags(mol, asMols=True)
                mol = max(frags, key=lambda m: m.GetNumAtoms())
                # Need to neutralize function group after remove charge
            if remove_stereo:
                Chem.RemoveStereochemistry(mol)
            
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol, canonical=True)
        
        except Exception as e:
            print(f"Error in standardizing {smiles}: {str(e)}")
            return mol 
        

    def smile_standardizer(self,desalt=True,remove_stereo=True):
        tqdm.pandas()
        self.df['standardized_smiles'] = self.df[self.smiles_column].progress_apply(lambda smiles: self.standardize_smiles(smiles,desalt=desalt,remove_stereo=remove_stereo), 
                                                                                           desc="Standardizing SMILES")

    def get_properties(self):
        self.df['mol'] = self.df['smiles'].apply(Chem.MolFromSmiles)
        self.df['MolWt'] = self.df['mol'].apply(Descriptors.MolWt)
        self.df['LogP'] = self.df['mol'].apply(Descriptors.MolLogP)
        self.df['NumHDonors'] = self.df['mol'].apply(Descriptors.NumHDonors)
        self.df['NumHAcceptors'] = self.df['mol'].apply(Descriptors.NumHAcceptors)
        self.df['TPSA'] = self.df['mol'].apply(Descriptors.TPSA)
        self.df['NumRotatableBonds'] = self.df['mol'].apply(Descriptors.NumRotatableBonds)
        self.df['RingCount'] = self.df['mol'].apply(Descriptors.RingCount)
        self.df['HeavyAtomCount'] = self.df['mol'].apply(Descriptors.HeavyAtomCount)
        self.df['FractionCSP3'] = self.df['mol'].apply(Descriptors.FractionCSP3)
        self.df['FormalCharge'] = self.df['mol'].apply(Chem.GetFormalCharge)
    
    def scaling(self):
        pass

    def pca_decomposition(self):
        pass

    def cluster(self):
        pass

    def sampling(self):
        pass
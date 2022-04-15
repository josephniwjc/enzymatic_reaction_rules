#! /usr/bin/env python

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
import pickle
import rdkit_utils


class PickXGBClassifier:
    """XGBoost model to predict feasibility of novel enzymatic reactions enumerated by Pickaxe"""

    def __init__(self, model_path, rules_path):
        """
        Parameters
        ----------
        model_path : str
            path to pickled enzymatic reaction feasibility classifier
        rules_path : str
            path to JN1224min ruleset
        """

        self.model = pickle.load(open(model_path, 'rb'))
        self.rules_df = pd.read_csv(rules_path, sep='\t', index_col=0)
        self._bondchange_dict = {}
        self.bondchange_featurization = lambda s: self._bondchange_lambda(s)
        self.compound_featurization = lambda s: self._compound_lambda(s)

    def predict_feasibility(self, reactant, product, rule, cutoff=0.5, return_proba=False):
        """
        Return feasibility or feasibility score of novel enzymatic reactions

        Parameters
        ----------
        reactant : str
            reactant SMILES
        product : str
            product SMILES
        rule : str
            rule id (for example 'rule0002')
        cutoff : float
            default 0.5, feasibility score above cutoff will be considered as feasible
        return_proba : bool
            default False, return feasibility score instead of feasibility if True
        """

        reaction_bits = np.hstack([self.bondchange_featurization(rule), self.compound_featurization(reactant),
                                   self.compound_featurization(product)])
        feasibility_score = self.model.predict_proba(reaction_bits.reshape(1, 5120))[0][1]

        if return_proba is True:
            return feasibility_score
        else:
            if feasibility_score >= cutoff:
                return True
            else:
                return False

    def _bondchange_lambda(self, rule):
        """Featurize bond change patterns"""

        try:
            rxn_array = self._bondchange_dict[rule]
        except KeyError:

            # extract bond change patterns
            lhs_smarts = self.rules_df.loc[rule, 'SMARTS'].split('>>')[0]
            rhs_smarts = self.rules_df.loc[rule, 'SMARTS'].split('>>')[1]
            lhs_any = rdkit_utils.get_smarts(lhs_smarts)[0][
                self.rules_df.loc[rule, 'Reactants'].split(';').index('Any')]
            rhs_any = rdkit_utils.get_smarts(rhs_smarts)[0][self.rules_df.loc[rule, 'Products'].split(';').index('Any')]

            # ECFP4
            lhs_ecfp4 = np.array([int(fp) for fp in DataStructs.cDataStructs.BitVectToText(
                Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(lhs_any), 4, nBits=256))])
            rhs_ecfp4 = np.array([int(fp) for fp in DataStructs.cDataStructs.BitVectToText(
                Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(rhs_any), 4, nBits=256))])
            lhs_ap = np.zeros(256)
            smarts_nonzero_elements = Chem.rdMolDescriptors.GetHashedAtomPairFingerprint(Chem.MolFromSmiles(lhs_any),
                                                                                         nBits=256).GetNonzeroElements()

            # Atom Pair
            for k, v in smarts_nonzero_elements.items():
                lhs_ap[k] = v
            rhs_ap = np.zeros(256)
            smarts_nonzero_elements = Chem.rdMolDescriptors.GetHashedAtomPairFingerprint(Chem.MolFromSmiles(rhs_any),
                                                                                         nBits=256).GetNonzeroElements()
            for k, v in smarts_nonzero_elements.items():
                rhs_ap[k] = v
            rxn_array = np.hstack([lhs_ecfp4, rhs_ecfp4, lhs_ap, rhs_ap])

            # store bond change fingerprint
            self._bondchange_dict[rule] = rxn_array

        return rxn_array

    def _compound_lambda(self, smiles):
        """Featurize reactant or product"""

        smiles_fp_array = np.zeros(2 * 1024, dtype=float)

        # ECFP4
        smiles_fp_array[0:1024] = [int(fp) for fp in DataStructs.cDataStructs.BitVectToText(
            Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 4, nBits=1024))]

        # Atom Pair
        smiles_nonzero_elements = Chem.rdMolDescriptors.GetHashedAtomPairFingerprint(Chem.MolFromSmiles(smiles),
                                                                                     nBits=1024).GetNonzeroElements()
        for k, v in smiles_nonzero_elements.items():
            smiles_fp_array[k + 1024] = v

        return list(smiles_fp_array)


if __name__ == '__main__':

    input_model_path = 'PickXGB_model.dat'
    input_rules_path = '../JN1224MIN/JN1224MIN_rules.tsv'
    PX = PickXGBClassifier(input_model_path, input_rules_path)

    print('Pyruvate to Alanine')
    print(PX.predict_feasibility('CC(=O)C(=O)O', 'CC(N)C(=O)O', 'rule0031', return_proba=True))
    print(PX.predict_feasibility('CC(=O)C(=O)O', 'CC(N)C(=O)O', 'rule0031'))
    print(PX.predict_feasibility('CC(=O)C(=O)O', 'CC(N)C(=O)O', 'rule0031', cutoff=0.970))
    print()

    print('Alanine to Betaalanine')
    print(PX.predict_feasibility('CC(N)C(=O)O', 'NCCC(=O)O', 'rule0028', return_proba=True))
    print(PX.predict_feasibility('CC(N)C(=O)O', 'NCCC(=O)O', 'rule0028'))
    print(PX.predict_feasibility('CC(N)C(=O)O', 'NCCC(=O)O', 'rule0028', cutoff=0.998))
    print()

    print('Betaalanine to Acrylate')
    print(PX.predict_feasibility('NCCC(=O)O', 'C=CC(=O)O', 'rule0244', return_proba=True))
    print(PX.predict_feasibility('NCCC(=O)O', 'C=CC(=O)O', 'rule0244'))
    print(PX.predict_feasibility('NCCC(=O)O', 'C=CC(=O)O', 'rule0244', cutoff=0.944))
    print()

    print('Acrylate to Propionate')
    print(PX.predict_feasibility('C=CC(=O)O', 'CCC(=O)O', 'rule0027', return_proba=True))
    print(PX.predict_feasibility('C=CC(=O)O', 'CCC(=O)O', 'rule0027'))
    print(PX.predict_feasibility('C=CC(=O)O', 'CCC(=O)O', 'rule0027', cutoff=0.997))
    print()

    print('Malate to 3-hydroxypropanoate')
    print(PX.predict_feasibility('O=C(O)CC(O)C(=O)O', 'O=C(O)CCO', 'rule0024', return_proba=True))
    print(PX.predict_feasibility('O=C(O)CC(O)C(=O)O', 'O=C(O)CCO', 'rule0024'))
    print(PX.predict_feasibility('O=C(O)CC(O)C(=O)O', 'O=C(O)CCO', 'rule0024', cutoff=0.180))

from rdkit import Chem
from rdkit.Chem import MolStandardize
from tqdm import tqdm
import re

class Preprocess:

    def __init__(self):
        self.nrm = MolStandardize.normalize.Normalizer()
        self.lfc = MolStandardize.fragment.LargestFragmentChooser()
        self.uc = MolStandardize.charge.Uncharger()
        self.pattern = re.compile("(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])")

    def clean(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = self.nrm.normalize(mol)
            mol = self.lfc.choose(mol)
            mol = self.uc.uncharge(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
        else:
            return None

    def tokenize(self, smi):
        return [token for token in self.pattern.findall(smi)]

if __name__ == "__main__":

    with open('data/canonical_smiles.smi', 'r') as file:
        smiles = [line.rstrip() for line in file]

    print("Initial number of sequences %i" % len(smiles))
    p = Preprocess()

    # Normalization, uncharging, removing chirality and light fragments
    nn_smi = [p.clean(smile) for smile in tqdm(smiles)]
    unn_smi = list(set([smile for smile in nn_smi if smile]))

    # Limit sequence length 34-74
    cl_smi = []
    for smile in unn_smi:
        if 34 <= len(p.tokenize(smile)) <= 74:
            cl_smi.append(smile)

    print("Number of sequences after cleaning %i" % len(cl_smi))

    with open('data/cleaned_smiles.smi', 'w') as file:
        for line in cl_smi:
            file.write(line + '\n')





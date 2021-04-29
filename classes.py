import keras
from keras.utils import to_categorical
import numpy as np
import re


class Tokenizer:

    def __init__(self):
        self.pattern = re.compile(
            '(A[ls]?|K|Br?|Cl?|Na?|O|S[ie]?|Te?|H|Zn?|P|F|I|b|c|n|o|se?|te?|p|\(|\)|\[|\]|\.|=|#|-|\+|\\|\/|:|~|@|\?|>|\*|\$|\%|[0-9]|G|A|E)'
        )

    def tokenize(self, sequence):
        return self.pattern.findall(sequence)

    def fit_on_texts(self, texts):
        if not isinstance(texts, list):
            texts = [texts]
        sequences = [self.tokenize(text) for text in texts]
        chars = list(set([symb for sequence in sequences for symb in sequence]))

        self.word_index = {val: i for i, val in enumerate(chars)}
        self.vocab_len = len(self.word_index)

    def texts_to_vector(self, texts):
        if not isinstance(texts, list):
            texts = [texts]
        try:
            return [[self.word_index[char] for char in self.tokenize(text)] for text in texts]
        except AttributeError as error:
            print('Object is not fitted yet. Call `fit_on_texts` instead.\n', error)


class DataGenerator(keras.utils.Sequence):

    def __init__(self, X, y, batch_size, num_classes):
        self.X = X
        self.y = y
        self.batch_size = batch_size
        self.num_classes = num_classes

    def __len__(self):
        return (np.ceil(len(self.X) / float(self.batch_size))).astype(np.int)

    def __getitem__(self, idx):
        batch_x = self.X[idx * self.batch_size: (idx + 1) * self.batch_size]
        batch_y = self.y[idx * self.batch_size: (idx + 1) * self.batch_size]

        return np.asarray([to_categorical(x, self.num_classes) for x in batch_x]), np.asarray(
            [to_categorical(x, self.num_classes) for x in batch_y])
